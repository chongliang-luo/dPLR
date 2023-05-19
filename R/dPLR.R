


#' @useDynLib stratPLR
#' @import stats
#' @import Rcpp
#' @title stratPLR: stratified Proportional Likelihood Ratio Model
#'
#' @description Fit a stratified PLR to integrate data across clinical sites.
#' The stratPLR likelihood does not involve the nuisance baseline distribution
#' functions, and hence allow the baseline distributions to be different from
#' sites. See Luo, et al. 2023 for details.
#' @usage stratPLR(mydata)
#' @author Chongliang Luo
#'
#' @param mydata Data frame of the data from all sites, with the first three columns as site id, outcome,
#'                 offset (e.g. observation time) and the rest columns as covariates.
#'
#' @return {an \code{optim()} output, optimizing the stratified PLR likelihood function.}
#'
#' @references Chongliang Luo, et al. "Distributed Proportional Likelihood Ratio Model with Application to Data Integration across Clinical Sites." AoAS (2023).
#'             Xiaodong Luo and Wei Yann Tsai. "A proportional likelihood ratio model." Biometrika (2012).
#' @export
stratPLR <- function(mydata){
  id.site.unique <- unique(mydata$id.site)
  px = ncol(mydata) - 3
  fn <- function(b){
    logL <- 0
    N2 <- 0
    for(ik in seq_along(id.site.unique)){
      sk <- id.site.unique[ik]
      di <- as.matrix(mydata[mydata$id.site==sk,-1])
      di <- di[order(di[,1]),]
      nk2 <- nrow(di) * (nrow(di)-1) / 2
      N2 <- N2 + nk2
      logL <- logL + nk2 * PropLikRatio(beta=b, y=di[,1], X=di[,-c(1,2)], z=di[,2], grad=0)$logL
      # logL <- logL + nrow(di) * PropLikRatio(beta=b, y=di[,1], X=di[,-c(1,2)], z=di[,2], grad=0)$logL
    }
    return(logL / N2)
  }

  fit_stratPLR <- optim(par = rep(0,px), fn = fn, hessian = T)
  return(fit_stratPLR)
}



#' @useDynLib dPLR
#' @import stats
#' @import Rcpp
#' @title dPLR: Distributed Proportional Likelihood Ratio Model
#'
#' @description Fit a privacy-preserving and communication-efficient dPLR model to integrate data across Clinical Sites. See Luo, et al. 2023 for details.
#'              This function shows the estimation accuracy of the method: dPLR is closer to the stratified PLR (pooled analysis) compared to the meta-estimator.
#' @usage dPLR(mydata, id.local=1, init_est = 'meta',  output.dPLR1 = F, control=list(maxit=100), verbose = F)
#' @author Chongliang Luo
#'
#' @param mydata Data frame of the data from all sites, with the first three columns as site id, outcome,
#'                 offset (e.g. observation time) and the rest columns as covariates.
#' @param id.local local site id
#' @param init_est Character or numeric vector, default 'meta' which uses PLR meta estimate as initial beta_bar.
#' @param output.dPLR1 Logical, if TRUE the surrogate likelihood is the 1st order approximation and PLR1 is also fitted. Default FALSE.
#' @param control control options passed to \code{optim()}, e.g. list(maxit=100).
#' @param verbose Logical - show running details?
#'
#' @return List with components:
#' @return \item{beta_bar}{initial coef estimate, by local site estimate or meta estimate}
#' @return \item{res_plr}{PLR fitted at all sites}
#' @return \item{b_meta}{meta PLR estimate}
#' @return \item{b_local}{local site PLR estimate}
#' @return \item{beta_tilde}{coef estimate by optimizing the distPLR surrogate likelihood function (2nd order approximation)}
#' @return \item{var_tilde}{var estimate of beta_tilde calculated from hessian}
#' @return \item{beta_tilde1}{coef estimate by optimizing the distPLR surrogate likelihood function (1st order approximation)}
#' @return \item{var_tilde1}{var estimate of beta_tilde1 calculated from hessian}
#' @return \item{sol}{an \code{optim()} output, optimizing distPLR surrogate likelihood function (2nd order approximation)}
#' @return \item{sol1}{an \code{optim()} output, optimizing distPLR surrogate likelihood function (1st order approximation)}
#' @return \item{grad_LN_L1_beta_bar}{sum gradients from collaborative sites, used in constructing the distPLR surrogate likelihood}
#' @return \item{hess_LN_L1_beta_bar}{sum hessians from collaborative sites, used in constructing the distPLR surrogate likelihood}
#' @details The surrogate pairwise log-likelihood is the logL of local site plus
#'          the 1st and 2nd order approximation of the logL of remote sites
#'          Ltilde(beta) = L1(beta) + sum(beta * grad(LN(beta_bar) - L1(beta_bar)))+ t(beta) %*% hessian(LN(beta_bar) - L1(beta_bar)) %*% beta / 2.
#'          It's optimized by \code{optim()}. The log-L, gradients and hessian are written in Rcpp.
#'
#' @references Chongliang Luo, et al. "Distributed Proportional Likelihood Ratio Model with Application to Data Integration across Clinical Sites." AoAS (2023).
#'             Xiaodong Luo and Wei Yann Tsai. "A proportional likelihood ratio model." Biometrika (2012).
#'              Michael I. Jordan, Jason D. Lee, and Yun Yang. "Communication-efficient distributed statistical inference." JASA (2018).
#' @examples
#' data(sim_data_dPLR)
#' ## distPLR
#' fit_dist <- dPLR(mydata=sim_data_dPLR, id.local = 1, init_est = 'meta')
#' b_plr_dist <- c(fit_dist$beta_tilde)
#' b_plr_meta <- fit_dist$b_meta
#'
#' ## pooling data together: stratified PLR (as gold-standard)
#' fit_strat = stratPLR(mydata=sim_data_dPLR)
#' b_plr_strat <- fit_strat$par
#'
#' ## dPLR is closer to gold-standard compared to the meta-estimator
#' round(cbind(b_plr_strat, b_plr_meta, b_plr_dist), 3)
#' @export
dPLR <- function(mydata,
                  id.local = 1,
                  init_est = 'local',
                  output.dPLR1 = FALSE,
                  control = list(maxit=300),
                  verbose = F){
  mydata$id.site <- as.character(mydata$id.site)
  id.site <- mydata$id.site
  id.site.unique <- names(sort(table(id.site), decreasing = T))

  px <- ncol(mydata) - 3
  if(!(id.local %in% id.site.unique)){
    cat('id.local not exist, use site', id.site.unique[1], 'as local')
    id.local <- id.site.unique[1]
  }

  ## local (single PLR) and meta (average PLR)
  K <- length(id.site.unique)
  nn <- table(id.site)
  N <- length(id.site)
  sum_K_beta <- rep(0, px)
  sum_K_wt <- matrix(0, px, px)  # rep(0, px)
  sum_K_beta_w <- rep(0, px)

  res_plr <- list()
  for(ik in id.site.unique){
    di <- mydata[mydata$id.site %in% c(ik), -1]
    di <- as.matrix(di[order(di[,1]),])
    fit_i <- PropLikRatio(y = di[,1], X = di[,-c(1,2)], z=di[,2], hessian=T)
    if(ik == id.local) b_local <- fit_i$par

    ## hessian as weight
    # wt <- fit_i$hessian
    wt <- diag(rep(nrow(di), px))
    sum_K_wt <- sum_K_wt + wt
    sum_K_beta_w <- sum_K_beta_w + wt %*% fit_i$par
    res_plr[[ik]] <- list(beta=fit_i$par, hessian=fit_i$hessian)
  }

  b_meta <-  solve(sum_K_wt, sum_K_beta_w)
  res_plr$meta <- b_meta


  # init est: beta_bar using local or meta
  if(is.numeric(init_est)){
    beta_bar <- init_est
  } else {
    if(init_est=='local'){
      if(verbose==T) cat('use the first site as init!')
      beta_bar <-  b_local
    } else if(init_est=='meta'){
      if(verbose==T) cat('use the meta PLR est as init!')
      beta_bar <- b_meta
    }
  }

  # logL and its gradient are written in rcpp
  logL_N_D1_beta_bar <- 0
  logL_N_D2_beta_bar <- 0
  for(ik in id.site.unique){  # setdiff(id.site, id.local)
    di <- mydata[mydata$id.site %in% c(ik), -1]
    di <- as.matrix(di[order(di[,1]),])
    nk <- nrow(di)
    gg <- PropLikRatio(beta=beta_bar, y = di[,1], X = di[,-c(1,2)], z=di[,2])
    logL_N_D1_beta_bar <- logL_N_D1_beta_bar + gg$logL_D1 * nk *(nk-1)/2
    logL_N_D2_beta_bar <- logL_N_D2_beta_bar + gg$logL_D2 * nk *(nk-1)/2
  }

  # surr estimates
  ik = id.local
  di <- mydata[id.site %in% c(ik), -1]
  di <- as.matrix(di[order(di[,1]),])
  nk <- nrow(di)

  logL <- function(b) PropLikRatio(beta=b, y = di[,1], X = di[,-c(1,2)], z=di[,2], grad=0)$logL
  # local derivatives
  tmp <- PropLikRatio(beta = beta_bar, y = di[,1], X = di[,-c(1,2)], z=di[,2])
  logL_D1_beta_bar <- tmp$logL_D1
  logL_D2_beta_bar <- tmp$logL_D2

  grad_LN_L1_beta_bar <- logL_N_D1_beta_bar / sum(nn*(nn-1)/2) - logL_D1_beta_bar
  hess_LN_L1_beta_bar <- logL_N_D2_beta_bar / sum(nn*(nn-1)/2) - logL_D2_beta_bar
  hess_LN_L1_beta_bar <- matrix(hess_LN_L1_beta_bar, px, px)

  # surrogate log-L and its gradient
  Ltilde1 <- function(b) logL(b) + sum(b * grad_LN_L1_beta_bar)
  Ltilde <- function(b) logL(b) + sum(b * grad_LN_L1_beta_bar) + 1/2 * t(c(b-beta_bar)) %*% hess_LN_L1_beta_bar %*% c(b-beta_bar)

  # Ltilde_gradient1 <- function(b) logL_gradient(b) + grad_LN_L1_beta_bar
  # Ltilde_gradient <- function(b) logL_gradient(b) + grad_LN_L1_beta_bar + hess_LN_L1_beta_bar %*% (b - beta_bar)

  if(output.dPLR1){
    sol1 <- optim(par = beta_bar, fn = Ltilde1, hessian = T, control = control)
    beta_tilde1 <- sol1$par
    var_tilde1 <- tryCatch(diag(solve(sol1$hessian))/N, error=function(e) NA)
  } else {
    sol1 <- list()
    beta_tilde1 <- NULL
    var_tilde1 <- NULL
  }

  sol <- optim(par = beta_bar, fn = Ltilde, hessian = T, control = control)
  if(verbose) cat('Ltilde at site ', ik, '...')

  beta_tilde <- sol$par

  # calculate var from inv hessian
  var_tilde <- tryCatch(diag(solve(sol$hessian))/N, error=function(e) NA)

  return(list(beta_bar = beta_bar,
              res_plr = res_plr,
              b_meta = c(b_meta),
              b_local = c(b_local),
              beta_tilde = c(beta_tilde),
              var_tilde  = var_tilde,
              beta_tilde1 = beta_tilde1,
              var_tilde1 = var_tilde1,
              sol = sol,
              sol1 = sol1,
              grad_LN_L1_beta_bar = grad_LN_L1_beta_bar,
              hess_LN_L1_beta_bar = hess_LN_L1_beta_bar ))
}




## PLR likelihood function. If beta is provided then compute the 1st and 2nd order
## derivatives, if not then optimize the PLR likelihood.
PropLikRatio <- function(beta=NULL, y, X, z=0, grad=NULL, hessian=F){
  # (-1/n) pairwise logL and its 1st, 2nd gradients,   Luo2012Biometrika
  # y: count data, need to be ordered
  # y.order <- order(y)
  X <- as.matrix(X)
  px <- ncol(X)
  index = c(0, cumsum(table(y)))
  if(is.null(beta)){
    if(is.null(grad)) grad <- 0
    mylogL <- function(b) PairwiseLik(yvec = y, z=z, X = as.matrix(X), index = index, beta = b, grad=grad)$logL
    res <- optim(par = rep(0, px), fn = mylogL, hessian=hessian)
  } else{
    if(is.null(grad)) grad <- 1
    res <- PairwiseLik(yvec = y, X = as.matrix(X), z=z, index = index, beta = beta, grad=grad)
  }
  return(res)
}


PLR.var <- function(bhat, mydata, option=1){
  id.site.unique <- unique(mydata$id.site)
  px <- length(bhat)
  nn = nn2 = nn3 = rep(NA, length(id.site.unique))

  Sigma_12_all <- list()
  S1_sum <- S2_sum <- matrix(0, px, px)
  for(ik in seq_along(id.site.unique)){
    sk <- id.site.unique[ik]
    di <- as.matrix(mydata[mydata$id.site==sk,-1])
    di <- di[order(di[,1]),]
    nn[ik] <- nrow(di)
    nn2[ik] <- nn[ik] * (nn[ik]-1) / 2
    # nn3[ik] <- nn[ik] * (nn[ik]-1) * (nn[ik]-1)/4
    # nn3[ik] <- nn[ik] * nn[ik] * (nn[ik]-1)/4

    if(option==0){
      S12 <- PairwiseLikDeriv(yvec = di[,1], z=di[,2], X = di[,-c(1,2)], beta = bhat)
      nn3[ik] <- nn[ik] * (nn[ik]-1)^2/4
    }else if(option==1){
      S12 <- PairwiseLikDeriv1(yvec = di[,1], z=di[,2], X = di[,-c(1,2)], beta = bhat)
      nn3[ik] <- nn[ik]^2 * (nn[ik]-1)/4
    }else if(option==2){
      S12 <- PairwiseLikDeriv2(yvec = di[,1], z=di[,2], X = di[,-c(1,2)], beta = bhat)
      nn3[ik] <- nn[ik] * (nn[ik]-1) * (nn[ik]-2) /4
    }

    Sigma_12_all[[as.character(sk)]] <- S12
    S1_sum <- S1_sum + S12$Sigma_1*nn2[ik]
    S2_sum <- S2_sum + S12$Sigma_2*nn3[ik]
  }
  S1_sum <- S1_sum / sum(nn2)
  S2_sum <- S2_sum / sum(nn3)
  V <- solve(S1_sum, S2_sum) %*% solve(S1_sum) / sum(nn)

  return(list(nn=nn, nn2=nn2, nn3=nn3, Sigma_12_all=Sigma_12_all,
              S1_sum=S1_sum, S2_sum=S2_sum, V=V))

}


## Luo2012Biometrika, pg 7, estimate PLR baseline g() given bhat
PLR.baseline <- function(bhat, mydata, eps=1e-5, max.iter=200){
  y <- mydata$y
  if(max(y) > 50) warning('max(y) > 50, numeric issue of eyXb may exist!')
  X <- as.matrix(mydata[,-c(1,2)])  # assume the first 2 cols are y and offset
  # z <- mydata[,2]
  z <- mydata$z
  n <- length(y)
  yu <- sort(unique(y))  # uniq value of y
  ny <- as.numeric(table(y))
  eyXb <-  exp(c(X%*%bhat + z) %*% t(yu) )
  eyXb[eyXb > 1e100] <- 1e100

  g0 <- ny / n
  g1 <- g0
  g.trace <- g0
  dif <- 1
  iter <- 1
  while(dif > eps & iter < max.iter){
    g0 <- g1
    A0 <-  c(eyXb %*% g0)
    g1 <- ny / apply(eyXb / A0, 2, sum)
    g1 <- g1 / sum(g1)
    dif <- max(abs(g1 - g0))
    g.trace <- rbind(g.trace, g1)
    iter <- iter + 1
    # cat('.')
  }

  conv <- TRUE
  if(iter == max.iter){
    conv <- FALSE
    cat('not converged yet, max(diff) =', dif)
  }
  return(list(g=data.frame(y=yu, prob=g1),
              g.trace=g.trace,
              iter=iter,
              conv=conv))
  # n <- 10000
  # X <- matrix(rnorm(n*2), n, 2)
  # bhat <- c(1,1)
  # y <- rpois(n, exp(X%*%bhat))
  # tmp = PLR.baseline(bhat=bhat, mydata=data.frame(y=y,0,X=X))
  # cbind(round(dpois(tmp$g$y, 1),4), round(tmp$g,4))
}


ztnb.dispersion <- function(y, X, z, fit){
  # z = offset
  if(any(y==0)){  # use (truncated) count part to calculate dispersion
    idx <- y!=0
    X=X[idx,]
    y=y[idx]
    z=z[idx]
    bhat = fit$coef$count
    yfit <- fit$fitted.values[idx]
  }else{
    bhat = fit$coef
    yfit <- fit$fitted.values
  }

  n <- length(y)
  px <- ncol(X)
  eXb <- c(exp(X%*% bhat[-1] + bhat[1] + z))
  r <- fit$theta
  p1 <- 1-(r/(r+eXb))^r
  yhat = eXb / p1
  vyhat <- (eXb + (1+1/r)*eXb^2)/p1 - (eXb/p1)^2
  resid <- (y - yhat) / sqrt(vyhat)
  phi.hat = sum(resid^2) / (n-px-1)
  resid <- (y - yfit) / sqrt(vyhat)
  phi.fit = sum(resid^2) / (n-px-1)
  return(c(phi.hat, phi.fit))
}

mspe <- function(y, X, z, bhat, theta, ghat, model=c('plr', 'ztnb', 'hurdle-nb')){
  n <- length(y)
  px <- ncol(X)
  if(model=='hurdle-nb'){
    bhat.count <- bhat[1:(px+1)]
    bhat.zero <- bhat[-c(1:(px+1))]
    eXb.count <- c(exp(X%*% bhat.count[-1] + bhat.count[1] + z))
    eXb.zero <- c(exp(X%*% bhat.zero[-1] + bhat.zero[1]))
    yu <- c(0:max(y))
    ydiff2 <- (outer(rep(1,n), yu) - y)^2
    prob <- matrix(NA, length(y), length(yu))
    prob[, 1] <- 1/(1+eXb.zero)
    for(iy in 1:max(y)) prob[, iy+1] <- 1/(1+1/eXb.zero) * dztnbinom(iy, eXb.count, theta)
    # prob <- prob / apply(prob,1,sum)
  } else if(model=='ztnb'){
    # bhat <- res$ztp$meta[1:(px+1)]
    eXb <- c(exp(X%*% bhat[-1] + bhat[1] + z))
    yu <- c(1:max(y))
    ydiff2 <- (outer(rep(1,n), yu) - y)^2
    prob <- matrix(NA, length(y), length(yu))
    for(iy in yu) prob[, iy] <- dztnbinom(iy, eXb, theta)
    # prob <- prob / apply(prob,1,sum)
  } else if(model=='plr'){
    # bhat <- res$dist$init_plr$beta_tilde
    eXb <- c(exp(X%*% bhat + z))
    yu <- ghat$y
    ydiff2 <- (outer(rep(1,n), yu) - y)^2
    # prob = outer(eyXb, ghat$prob)
    prob <- matrix(NA, length(y), length(yu))
    for(iy in seq_along(ghat$y)) prob[, iy] <- eXb^(ghat$y[iy]) * ghat$prob[iy]
    # prob <- prob / apply(prob,1,sum)
  }

  # mean square predictive error: sum_k (y-k)^2 * Pr(y=k)
  mspe0 <- list()
  pp <- prob / apply(prob,1,sum)
  mspe0[['y']] <- apply(ydiff2 * pp, 1, sum)
  mspe0[['ymae']] <- apply(sqrt(ydiff2) * pp, 1, sum)
  ## truncated MSPE...
  pp <- prob[,yu<=10] / apply(prob[,yu<=10],1,sum)
  mspe0[['y10']] <- apply(ydiff2[y<=10, yu<=10] * pp[y<=10,], 1, sum)
  mspe0[['y10mae']] <- apply(sqrt(ydiff2[y<=10, yu<=10]) * pp[y<=10,], 1, sum)
  pp <- prob[,yu<=5] / apply(prob[,yu<=5],1,sum)
  mspe0[['y5']] <- apply(ydiff2[y<=5, yu<=5] * pp[y<=5,], 1, sum)
  mspe0[['y5mae']] <- apply(sqrt(ydiff2[y<=5, yu<=5]) * pp[y<=5,], 1, sum)
  mspe0[['avg']] <- c(mean(mspe0[['y']]), mean(mspe0[['ymae']]),
                      mean(mspe0[['y10']]), mean(mspe0[['y10mae']]),
                      mean(mspe0[['y5']]), mean(mspe0[['y5mae']]))
  return(mspe0)
}
