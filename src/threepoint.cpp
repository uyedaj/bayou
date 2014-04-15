#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

using namespace arma; 
using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
SEXP C_threepoint(SEXP dat) 
  {
    Rcpp::List cache(dat);
    int root = Rcpp::as<int>(cache["root"]);
    int n = Rcpp::as<int>(cache["n"]);
    int N = Rcpp::as<int>(cache["N"]);
    int Nnode = N - n + 1;
    int ROOT = n ;
    arma::vec len = Rcpp::as<arma::vec >(cache["len"]);
    arma::vec diagm = Rcpp::as<arma::vec >(cache["diagMatrix"]);
    arma::vec P = Rcpp::as<arma::vec >(cache["P"]);
    std::vector<int> anc = Rcpp::as<std::vector<int> >(cache["anc"]);
    std::vector<int> des = Rcpp::as<std::vector<int> >(cache["des"]);
    int rowP = P.n_rows;
    int des_i, anc_i;
    arma::mat Pm;
    Pm.ones(rowP, 2);
    Pm.col(1) = P;
    arma::mat P_dm = Pm;
    P_dm.col(0) = P_dm.col(0) / diagm;
    P_dm.col(1) = P_dm.col(1) / diagm;
    arma::cube PP;
    PP.zeros(2, 2, n + Nnode);
    arma::vec zero, tmp11, logd;
    arma::mat P_dm_des, pp, tmpP1, tmp;
    tmpP1.zeros(2, n + Nnode);
    P_dm_des.zeros(1, 2);
    pp.zeros(2,2);
    zero = tmp11 = logd = zeros(n + Nnode);
    double el;
    //double tmp2;
    for (int i=0; i<N; ++i){
      el = len[i];
      des_i = des[i] - 1;
      anc_i = anc[i] - 1;
      if(des_i < n){
                if(el > 0){
                  logd[des_i] = log(el);
                  P_dm_des = P_dm.row(des_i);
                  pp = P_dm_des.t() * P_dm_des;
                  PP.subcube(0, 0, des_i, 1, 1, des_i) = pp / el;
                  tmpP1.col(des_i) = trans(P_dm.row(des_i) / el);
                  tmp11(des_i) = 1/el;  	  
                } else zero(anc_i) = zero(anc_i) + 1;
      } else {
        if((el <= 0 && zero(des_i) > 0) || zero(des_i) >1) break;
          logd(des_i) = logd(des_i) + log(1+el*tmp11(des_i));
          tmp =  tmpP1.col(des_i) * trans(tmpP1.col(des_i));
          //tmp2 = el/(1 + el * tmp11(des_i));
          tmp *= el/(1 + el * tmp11(des_i));
          PP.subcube(0, 0, des_i, 1, 1, des_i) += -1*tmp;
          tmpP1.col(des_i) *= (1 - 1/(1+1/el/tmp11(des_i)));     	
          tmp11(des_i) *= (1 - 1/(1+1/el/tmp11(des_i)));
      };
      logd(anc_i) += logd(des_i);
      PP.subcube(0, 0, anc_i, 1, 1, anc_i) += PP.subcube(0, 0, des_i, 1, 1, des_i);
      tmpP1.col(anc_i) += tmpP1.col(des_i);
      tmp11(anc_i) += tmp11(des_i);
    }
    logd(ROOT) += log(1+root*tmp11(ROOT)) + 2*sum(log(diagm));
    PP.subcube(0, 0, ROOT, 1, 1, ROOT) += -1*root/(1+root*tmp11(ROOT))*tmpP1.col(ROOT)*trans(tmpP1.col(ROOT));  
    arma::vec vec11, vecPP, vecP1;
    vec11 = PP.subcube(0,0,ROOT,0,0,ROOT);
    vecP1 = PP.subcube(0,1,ROOT,0,1,ROOT);
    vecPP = PP.subcube(1,1,ROOT,1,1,ROOT);
    return Rcpp::List::create(
      Rcpp::Named("P_dm", P_dm),
      Rcpp::Named("tmpP1", tmpP1),
      Rcpp::Named("vec11", vec11),
      Rcpp::Named("P1", vecP1),
      Rcpp::Named("PP", vecPP),
      Rcpp::Named("logd", logd(ROOT))
		);
  }
  
// [[Rcpp::export]]
SEXP C_transf_branch_lengths(SEXP dat, int model, NumericVector y, double alpha) 
  { 
    Rcpp::List cache(dat);
    int N = Rcpp::as<int>(cache["N"]);
    int n = Rcpp::as<int>(cache["n"]);
    LogicalVector externalEdge = Rcpp::as<LogicalVector>(cache["externalEdge"]);
    arma::vec times = Rcpp::as<arma::vec >(cache["times"]);
    arma::vec D = Rcpp::as<arma::vec >(cache["D"]);
    std::vector<int> banc = Rcpp::as<std::vector<int> >(cache["branches.anc2"]);
    std::vector<int> bdes = Rcpp::as<std::vector<int> >(cache["branches.des2"]);
    std::vector<int> anc = Rcpp::as<std::vector<int> >(cache["anc"]);
    std::vector<int> des = Rcpp::as<std::vector<int> >(cache["des"]);
    double Tmax = Rcpp::as<double>(cache["Tmax"]);
    //int anc_i, des_i, banc_i, bdes_i;
    int des_i, banc_i, bdes_i;
    double d1, d2, root;
    arma::vec len, distFromRoot, exp_times, diagMatrix;
    if(alpha==0){
      len = Rcpp::as<arma::vec >(cache["edge.length"]);
      root = 0;
      diagMatrix.ones(n);
    } else {
      if(model == 1){
      len.zeros(N);
      exp_times = exp(-2*alpha*times);
      distFromRoot = (1 - exp(-2 * alpha * (Tmax - times)));
      distFromRoot %= exp_times;
      for (int i=0; i<N; ++i){
      //anc_i = anc[i] - 1;
        des_i = des[i] - 1;
        banc_i = banc[i] - 1;
        bdes_i = bdes[i] - 1;
        d1 = distFromRoot(banc_i);
        if(externalEdge[i]) {
          d2 = exp(-2*alpha*D(des_i)) * (1-exp(-2*alpha*(Tmax-D(des_i))));
        } else {
          d2 = distFromRoot(bdes_i);
        };
        len(i) = d2 - d1;
      };
      root = distFromRoot.min();
      diagMatrix = exp(alpha*D);
      //edge.length = numeric(cache$N)
      //distFromRoot <-  exp(-2*alpha*cache$times)*(1 - exp(-2*alpha*(cache$Tmax-cache$times)))
      //#for (i in 1:cache$N) {
      //  #d1 = distFromRoot[which(names(cache$times) == cache$anc[i])]
      //#  d1 = distFromRoot[cache$branches.anc[[i]]]
      //#  if (cache$externalEdge[i]) {d2 = exp(-2*alpha*cache$D[cache$des[i]])*(1-exp(-2*alpha*(cache$Tmax-cache$D[cache$des[i]])))}
      //#  else #d2 = distFromRoot[which(names(cache$times) == cache$des[i])]
      //#/    {d2 = distFromRoot[cache$branches.des[[i]]]}
      //#  edge.length[i] = d2 - d1
      //#}
      //edge.length <- sapply(1:cache$N, function(i){
      //  d1 = distFromRoot[cache$branches.anc[[i]]];
      //  if (cache$externalEdge[i]) {d2 = exp(-2*alpha*cache$D[cache$des[i]])*(1-exp(-2*alpha*(cache$Tmax-cache$D[cache$des[i]])))} else d2 = distFromRoot[cache$branches.des[[i]]];
      //  d2 - d1
      //})
      //root.edge = min(distFromRoot)
      //diagMatrix = exp(alpha*cache$D)
      }
    }
    return Rcpp::List::create(
      Rcpp::Named("edge.length", len),
      Rcpp::Named("root.edge", root),
      Rcpp::Named("diagMatrix", diagMatrix)
    );
  }