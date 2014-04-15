#include <RcppArmadillo.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace arma; 
using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]  
SEXP C_weightmatrix(SEXP dat, SEXP parameters) 
  {
    Rcpp::List cache(dat);
    Rcpp::List pars(parameters);
    int k = Rcpp::as<int>(pars["k"]);
    int n = Rcpp::as<int>(cache["n"]);
    arma::vec theta = Rcpp::as<arma::vec >(pars["theta"]);
    arma::vec data = Rcpp::as<arma::vec >(cache["dat"]);
    if(k==0){
      arma::vec W(n);
      W.ones();
      arma::vec Etheta = W * theta;
      arma::vec residuals = data - Etheta;
      return Rcpp::List::create(
        Rcpp::Named("W", W),
        Rcpp::Named("E", Etheta),
        Rcpp::Named("resid", residuals)
      );
    } else {
    int N = Rcpp::as<int>(cache["N"]);
    int ntheta = Rcpp::as<int>(pars["ntheta"]);
    double alpha = Rcpp::as<double>(pars["alpha"]);
    bool ultrametric = Rcpp::as<bool>(cache["ultrametric"]);
    arma::ivec sb = Rcpp::as<arma::ivec >(pars["sb"]);
    arma::ivec t2 = Rcpp::as<arma::ivec >(pars["t2"]);
    arma::vec loc = Rcpp::as<arma::vec >(pars["loc"]);
    arma::vec len = Rcpp::as<arma::vec >(cache["edge.length"]);
    arma::vec nH = Rcpp::as<arma::vec >(cache["nH"]);
    arma::ivec anc = Rcpp::as<arma::ivec >(cache["anc"]);
    arma::ivec des = Rcpp::as<arma::ivec >(cache["des"]);
    arma::mat branchtrace = Rcpp::as<arma::mat >(cache["branchtrace"]);
    arma::vec tipheight;
    if(ultrametric > 0){
      tipheight = Rcpp::as<arma::vec >(cache["height"]);
    } else {
      tipheight = Rcpp::as<arma::vec >(cache["tipFromRoot"]);
    }
    arma::ivec start(N);
    start.ones();
    arma::uvec l_ord = stable_sort_index(loc);
    sb = sb(l_ord);
    t2 = t2(l_ord);
    loc = loc(l_ord);
    arma::mat bW(N, ntheta);
    bW.zeros();
    for(int i=N; i > 0; --i){
      //int i = 15;
      arma::vec mstart;
      arma::ivec mtheta(1);
      arma::vec mend(1);
      mstart.zeros(1);
      mend(0) = len(i-1);
      arma::uvec ind = find(sb==i);
      mtheta(0) = start(i-1);
      arma::vec x = loc(ind);
      if(ind.n_rows > 0){
        mstart.insert_rows(mstart.n_rows, loc(ind));
        mend.insert_rows(0, loc(ind));
        mtheta.insert_rows(mtheta.n_rows, t2(ind));
      }
        arma::uvec desb = find(anc == des(i-1));
        if(desb.n_rows >0){
          start(desb(0)) = (mtheta(mtheta.n_rows-1));
          start(desb(1)) = (mtheta(mtheta.n_rows-1));
        }
        arma::vec mdiff = mend - mstart;
        arma::vec W1 = alpha*(nH(i-1) + mstart - tipheight.max());
        arma::vec W2 = alpha * mdiff;
        for(int j=0; j < mstart.n_rows; ++j){
          if(any(W2 > 500)){
            if(W2(j) > 500){
              bW.submat(i-1, mtheta(j)-1, i-1, mtheta(j)-1) =  bW.submat(i-1, mtheta(j)-1, i-1, mtheta(j)-1) + exp(W1(j)+W2(j));
            } else{
              bW.submat(i-1, mtheta(j)-1, i-1, mtheta(j)-1) =  bW.submat(i-1, mtheta(j)-1, i-1, mtheta(j)-1) + exp(W1(j)) * expm1(W2(j));
            }
          } else {
            bW.submat(i-1, mtheta(j)-1, i-1, mtheta(j)-1) =  bW.submat(i-1, mtheta(j)-1, i-1, mtheta(j)-1) + exp(W1(j)) * expm1(W2(j));
          }
        }
      }
      arma::mat W(n, ntheta);
      arma::mat btrace;
      if(ultrametric==1){
        arma::mat btrace = branchtrace;
        W = btrace * bW;
        W.col(0) = W.col(0) + as_scalar(exp(-alpha * tipheight));
      } else {
        arma::mat btrace = branchtrace;
        for(int i=0; i < n; ++i){
          btrace.row(i) = btrace.row(i) * exp(alpha * (tipheight.max() - tipheight(i)));
        }
        W = btrace * bW;
        W.col(0) = W.col(0) + exp(-alpha * tipheight);
      }
    arma::vec Etheta = W * theta;
    arma::vec residuals = data - Etheta;
    return Rcpp::List::create(
      Rcpp::Named("bW", bW),
      Rcpp::Named("W", W),
      Rcpp::Named("E", Etheta),
      Rcpp::Named("resid", residuals)
      );
    }
  }
  
  
