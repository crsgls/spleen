#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat mspline(arma::colvec x, double tmin, double tmax, arma::colvec tint, int k){
  
  int lgtx = x.size();
  int lgtint = tint.size();
  
  arma::colvec t = arma::zeros<arma::colvec>(lgtint+2*k);
  arma::colvec tmin_vec = arma::zeros<arma::colvec>(k);
  arma::colvec tmax_vec = arma::zeros<arma::colvec>(k);
  arma::mat Mspl = arma::zeros<arma::mat>(lgtx, lgtint + k);
  
  
  for (int i=0; i < k; i++) {
    tmin_vec(i) = tmin; tmax_vec(i) = tmax;
  }
  t = join_cols(join_cols(tmin_vec, tint), tmax_vec);
  
  if (k == 1){
    for (int idx=0; idx<lgtx; idx++) {
      for (int i=0; i<k+lgtint; i++) {
        if((x(idx)>=t(i)) & (x(idx)<t(i+1))){
          Mspl(idx,i) = 1/(t(i+1)-t(i));
        }
      }
    }
  }
  
  if (k > 1){
    arma::mat mv = mspline(x,tmin,tmax,tint,k-1);
    for (int idx=0; idx<lgtx; idx++) {
      for (int i=0; i<k+lgtint; i++) {
        double mi = 0; double mii = 0;
        if((x(idx) < std::min(t(i+k-1),tmax)) & (x(idx) >= std::max(t(i),tmin))){
          mi = mv(idx,i-1);
        }
        if((x(idx) >= std::max(t(i+1),tmin)) & (x(idx) < std::min(t(i+k),tmax))){
          mii = mv(idx,i);
        }
        Mspl(idx,i) = k*((x(idx)-t(i))*mi + (t(i+k)-x(idx))*mii)/((k-1)*(t(i+k)-t(i)));
      }
    }
  }
  return(Mspl);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat ispline(arma::colvec x, double tmin, double tmax, arma::colvec tint, int k){
  
  int lgtx = x.size();
  int lgtint = tint.size();
  
  arma::colvec t = arma::zeros<arma::colvec>(lgtint+2*(k+1));
  arma::colvec tmin_vec = arma::zeros<arma::colvec>(k+1);
  arma::colvec tmax_vec = arma::zeros<arma::colvec>(k+1);
  arma::mat Ispl = arma::zeros<arma::mat>(lgtx, lgtint + k);
  
  for (int i=0; i<k+1; i++) {
    tmin_vec(i) = tmin; tmax_vec(i) = tmax;
  }
  
  t = join_cols(join_cols(tmin_vec, tint), tmax_vec);
  
  arma::mat mv = mspline(x, tmin, tmax, tint, k+1);
  for (int idx=0; idx<lgtx; idx++){
    
    int J = 0;
    for (int j = 0; j<lgtint+k+1; j++){
      if((x(idx)<t(j+1)) & (x(idx)>=t(j))){
        J = j;
      }
    }
    
    for (int i=1; i<k+lgtint+2; i++){
      if (i < J-k+1){
        Ispl(idx,i-1) = 1;
      }
      
      if((J-k+1 <= i) & (i <= J)){
        for (int l=i; l<J+1; l++){
          Ispl(idx,i-1) = Ispl(idx,i-1) + (t(l+k+1)-t(l))*mv(idx,l)/(k+1);
        }
      }
    }
  }
  
  return(Ispl);
}