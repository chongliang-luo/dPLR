#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <R.h>
#include <stdio.h>
#include <math.h>

using namespace arma;
using namespace Rcpp;
using namespace std;
// using namespace sugar;
// need to use older version of RcppArmadillo...
// install.packages("~/Dropbox/R/DistPoi/ODAP/RcppArmadillo_0.9.900.2.0.tar.gz", repos = NULL, type = "source")


// [[Rcpp::export]]
Rcpp::List PairwiseLik(arma::vec yvec, arma::colvec z, arma::mat X, arma::vec index, arma::colvec beta, int grad){
  // (-1/n) pairwise logL and its 1st, 2nd gradients,   Luo2012Biometrika
  // yvec is count outcome, need to be ordered
  // z is the offset term, e.g. the log(exposure time)
  // index = c(0, cumsum(table(y)))
  int i, j, ii, jj;
  int px = X.n_cols, n = X.n_rows, ngroup = index.n_elem - 1;
  int n2 = n * (n-1) / 2;
  arma::rowvec D1 = zeros<rowvec>(px);
  arma::mat D2(px, px, fill::zeros);

  double logL = 0;
  double eyXb, ydiff;
  arma::rowvec Xdiff;
  arma::colvec Xb = X*beta + z;

  for(i=0; i<(ngroup-1); i++){
    for(j=i+1; j<ngroup; j++){
      for(ii=index(i); ii<index(i+1); ii++){
        for(jj=index(j); jj<index(j+1); jj++){
          Xdiff = X.row(ii) - X.row(jj);
          ydiff = yvec(ii) - yvec(jj);
          eyXb = exp(-ydiff * (Xb(ii) - Xb(jj)));

          logL = logL - log(1+eyXb);
          if(grad==1){
            D1 = D1 + Xdiff * ydiff * (1 - 1 / (1+eyXb));
            // logL = 0;
            D2 = D2 - (Xdiff.t() * Xdiff) * pow(ydiff, 2) * eyXb / pow(1+eyXb, 2);
          }

        }
      }
    }
  }

  Rcpp::List out;
  out["logL"] = - logL / n2;
  if(grad==1){
    out["logL_D1"] = - D1 / n2;
    out["logL_D2"] = - D2 / n2;
  }
  return(out);
}




// [[Rcpp::export]]
Rcpp::List PairwiseLikDeriv(arma::vec yvec, arma::colvec z, arma::mat X, arma::colvec beta){
  // yvec is count outcome, need to be ordered
  // z is the offset term, e.g. the log(exposure time)
  // Var estimation: D2 = Sigma_1=E(dS(beta)/dbeta), D3 = Sigma_2 = Cov(S(beta))
  int i, ii, iii;
  int px = X.n_cols, n = X.n_rows;
  int n2 = n * (n-1) / 2 ;
  // int n3 = n * n * (n-1) / 2;
  int n3 = n * (n-1) * (n-1);
  // int n3 = n * (n-1) * (n-2);
  arma::mat D2(px, px, fill::zeros);
  arma::mat D3(px, px, fill::zeros);
  arma::mat phisq(px, px, fill::zeros);  // phi * phi^T
  double eyXb, ydiff;
  arma::rowvec Xdiff, Xdiff2, d1, d2, phi(px), phi12(px), phi13(px);
  arma::colvec Xb = X*beta + z;
  phi.fill(0);
  phi12.fill(0);
  phi13.fill(0);
  arma::cube phiall(n, n, px, fill::zeros); // phiall.fill(0);

  for(i=0; i<n; i++){
    for(ii=0; ii<n; ii++){
      if(ii>i){
        ydiff = yvec(i) - yvec(ii);
        Xdiff = X.row(i) - X.row(ii);
        eyXb = exp(-ydiff * (Xb(i) - Xb(ii)));
        D2 = D2 - (Xdiff.t() * Xdiff) * pow(ydiff, 2) * eyXb / pow(1+eyXb, 2);
        phi12 = (1-1/(1+eyXb))*ydiff*Xdiff;
        phiall.subcube(i, ii, 0, i, ii, px-1) = phi12;
        phi += phi12;
      }
      if(ii<i){
        phiall.subcube(i, ii, 0, i, ii, px-1) = phiall.subcube(ii, i, 0, ii, i, px-1);
      }
    }
  }
  // D2 is actually hessian from optim (at bhat), can be omitted...
  D2 = D2 / n2;
  // phi is the score (at bhat), ~0 can be omitted...
  phi = phi / n2;
  phisq = phi.t() * phi;
  // std::cout << "p4" << endl;

  for(i=0; i<n; i++){
    // for(ii=0; ii<(n-1); ii++){
      for(ii=0; ii<n; ii++){
      // for(iii=ii+1; iii<n; iii++){
        for(iii=0; iii<n; iii++){
          if(ii!=i & iii!=i){
            phi12 = phiall.subcube(i, ii, 0, i, ii, px-1);
            phi13 = phiall.subcube(i, iii, 0, i, iii, px-1);
            // D3 += 4*(phi12.t()*phi13 - phisq);
            D3 += 4*(phi12.t()*phi13);
          }
      }
    }
  }
  // D3 = D3 / n3;
  D3 = D3 / n3 - phisq;  //; //

  Rcpp::List out;
  out["Sigma_1"] = D2;
  out["Sigma_2"] = D3;
  out["Phi"] = phi;
  return(out);
}


// [[Rcpp::export]]
Rcpp::List PairwiseLikDeriv1(arma::vec yvec, arma::colvec z, arma::mat X, arma::colvec beta){
  // yvec is count outcome, need to be ordered
  // z is the offset term, e.g. the log(exposure time)
  // Var estimation: D2 = Sigma_1=E(dS(beta)/dbeta), D3 = Sigma_2 = Cov(S(beta))
  int i, ii; //, iii
  int px = X.n_cols, n = X.n_rows;
  int n2 = n * (n-1) / 2 ;
  // int n3 = n * n * (n-1) / 2;
  // int n3 = n * n * (n-1) / 2;
  arma::mat D2(px, px, fill::zeros);
  arma::mat D3(px, px, fill::zeros);
  arma::mat phisq(px, px, fill::zeros);  // phi * phi^T
  double eyXb, ydiff;
  arma::rowvec Xdiff, Xdiff2, d1, d2, phi(px), phi12(px), phi13(px);
  arma::colvec Xb = X*beta + z;
  phi.fill(0);
  phi12.fill(0);
  phi13.fill(0);
  arma::cube phiall(n, n, px, fill::zeros); // phiall.fill(0);

  for(i=0; i<n; i++){
    for(ii=0; ii<n; ii++){
      if(ii>i){
        ydiff = yvec(i) - yvec(ii);
        Xdiff = X.row(i) - X.row(ii);
        eyXb = exp(-ydiff * (Xb(i) - Xb(ii)));
        D2 = D2 - (Xdiff.t() * Xdiff) * pow(ydiff, 2) * eyXb / pow(1+eyXb, 2);
        phi12 = (1-1/(1+eyXb))*ydiff*Xdiff;
        phiall.subcube(i, ii, 0, i, ii, px-1) = phi12;
        phi += phi12;

      }
      if(ii<i){
        phiall.subcube(i, ii, 0, i, ii, px-1) = phiall.subcube(ii, i, 0, ii, i, px-1);
      }
    }
  }
  D2 = D2 / n2;
  phi = phi / n2;
  phisq = phi.t() * phi;
  // std::cout << "p4" << endl;

  for(i=0; i<n; i++){
    phi12.fill(0);
    for(ii=0; ii<n; ii++){
    // for(ii=0; ii<(n-1); ii++){
      // for(iii=ii+1; iii<n; iii++){
      if(ii != i){
          phi12 += phiall.subcube(i, ii, 0, i, ii, px-1);
      }
    }
    phi12 = phi12 / (n-1) - phi;
    D3 += 4*(phi12.t()*phi12);
          // phi13 = phiall.subcube(i, iii, 0, i, iii, px-1);
          // D3 += 4*(phi12.t()*phi13 - phisq);
      // }
    // }
  }
  D3 = D3 / (n-1);

  Rcpp::List out;
  out["Sigma_1"] = D2;
  out["Sigma_2"] = D3;
  out["Phi"] = phi;
  return(out);
}


// [[Rcpp::export]]
Rcpp::List PairwiseLikDeriv2(arma::vec yvec, arma::colvec z, arma::mat X, arma::colvec beta){
  // yvec is count outcome, need to be ordered
  // z is the offset term, e.g. the log(exposure time)
  // Var estimation: D2 = Sigma_1=E(dS(beta)/dbeta), D3 = Sigma_2 = Cov(S(beta))
  int i, ii, iii;
  int px = X.n_cols, n = X.n_rows;
  int n2 = n * (n-1) / 2 ;
  // int n3 = n * n * (n-1) / 2;
  // int n3 = n * (n-1) * (n-1);
  int n3 = n * (n-1) * (n-2);
  arma::mat D2(px, px, fill::zeros);
  arma::mat D3(px, px, fill::zeros);
  arma::mat phisq(px, px, fill::zeros);  // phi * phi^T
  double eyXb, ydiff;
  arma::rowvec Xdiff, Xdiff2, d1, d2, phi(px), phi12(px), phi13(px);
  arma::colvec Xb = X*beta + z;
  phi.fill(0);
  phi12.fill(0);
  phi13.fill(0);
  arma::cube phiall(n, n, px, fill::zeros); // phiall.fill(0);

  for(i=0; i<n; i++){
    for(ii=0; ii<n; ii++){
      if(ii>i){
        ydiff = yvec(i) - yvec(ii);
        Xdiff = X.row(i) - X.row(ii);
        eyXb = exp(-ydiff * (Xb(i) - Xb(ii)));
        D2 = D2 - (Xdiff.t() * Xdiff) * pow(ydiff, 2) * eyXb / pow(1+eyXb, 2);
        phi12 = (1-1/(1+eyXb))*ydiff*Xdiff;
        phiall.subcube(i, ii, 0, i, ii, px-1) = phi12;
        phi += phi12;
      }
      if(ii<i){
        phiall.subcube(i, ii, 0, i, ii, px-1) = phiall.subcube(ii, i, 0, ii, i, px-1);
      }
    }
  }
  D2 = D2 / n2;
  phi = phi / n2;
  phisq = phi.t() * phi;
  // std::cout << "p4" << endl;

  for(i=0; i<n; i++){
    // for(ii=0; ii<(n-1); ii++){
    for(ii=0; ii<n; ii++){
      // for(iii=ii+1; iii<n; iii++){
      for(iii=0; iii<n; iii++){
        if(ii!=i & iii!=i & ii!=iii){
          phi12 = phiall.subcube(i, ii, 0, i, ii, px-1);
          phi13 = phiall.subcube(i, iii, 0, i, iii, px-1);
          // D3 += 4*(phi12.t()*phi13 - phisq);
          D3 += 4*(phi12.t()*phi13);
        }
      }
    }
  }
  // D3 = D3 / n3;
  D3 = D3 / n3 - phisq;

  Rcpp::List out;
  out["Sigma_1"] = D2;
  out["Sigma_2"] = D3;
  out["Phi"] = phi;
  return(out);
}


// [[Rcpp::export]]
Rcpp::List PairwiseLikDeriv3(arma::vec yvec, arma::colvec z, arma::mat X, arma::colvec beta){
  // yvec is count outcome, need to be ordered
  // z is the offset term, e.g. the log(exposure time)
  // index = c(0, cumsum(table(y)))
  // Var estimation: D2 = Sigma_1=E(dS(beta)/dbeta), D3 = Sigma_2 = Cov(S(beta))
  int i, j, j2;
  int px = X.n_cols, n = X.n_rows;
  int n2 = n * (n-1) ;
  int n3 = n * (n-1) * (n-1) ;
  // arma::rowvec D1 = zeros<rowvec>(px);
  arma::mat D2(px, px, fill::zeros);
  arma::mat D3(px, px, fill::zeros);

  // double logL = 0;
  double eyXb, ydiff, eyXb2, ydiff2;
  arma::rowvec Xdiff, Xdiff2, d1, d2;
  arma::colvec Xb = X*beta + z;

  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      ydiff = yvec(i) - yvec(j);
      if(ydiff!=0){
        Xdiff = X.row(i) - X.row(j);
        eyXb = exp(-ydiff * (Xb(i) - Xb(j)));
        D2 = D2 - (Xdiff.t() * Xdiff) * pow(ydiff, 2) * eyXb / pow(1+eyXb, 2);

        for(j2=0; j2<n; j2++){
          ydiff2 = yvec(i) - yvec(j2);
          if(ydiff2!=0){
            Xdiff2 = X.row(i) - X.row(j2);
            eyXb2 = exp(-ydiff2 * (Xb(i) - Xb(j2)));
            d1 = Xdiff * ydiff * (1 - 1 / (1+eyXb));
            d2 = Xdiff2 * ydiff2 * (1 - 1 / (1+eyXb2));
            D3 = D3 - d1.t() * d2;
          }
        }
      }
    }
  }
  // std::cout << logL << D1 << D2 << endl;
  Rcpp::List out;
  //out["logL"] = - logL / n2;
  out["Sigma_1"] = - D2 / n2;
  out["Sigma_2"] = - D3 * 4 / n3;

  return(out);
}
