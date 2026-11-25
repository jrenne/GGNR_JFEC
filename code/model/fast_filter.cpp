// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include "pricing_helpers.h"
using namespace Rcpp;
static double filter_norm_cdf(double x) { return R::pnorm5(x,0.0,1.0,1,0); }
static double filter_norm_pdf(double x) { return R::dnorm4(x,0.0,1.0,0); }

struct Measurement {
  arma::mat G, delta, cross;
  arma::vec innovation;
};

// Same IEKF equations and pricing routine as the published R implementation.
// This implementation does not alter the model or differentiation rule.
// [[Rcpp::export]]
Rcpp::List fast_iekf_cpp(Rcpp::List obj, Rcpp::NumericMatrix observations,
 Rcpp::IntegerVector mats_n, Rcpp::IntegerVector mats_r,
 Rcpp::IntegerVector horizons, Rcpp::NumericVector initial_mean,
 Rcpp::NumericMatrix initial_covariance, bool iterated_filter,
 double trigger, int steps, Rcpp::Nullable<Rcpp::LogicalVector> branches,
 bool record) {
 using arma::mat; using arma::vec;
 mat Y=Rcpp::as<mat>(observations), Phi=Rcpp::as<mat>(obj["Phi"]),
   Sigma=Rcpp::as<mat>(obj["Sigma2_X"]), P=Rcpp::as<mat>(initial_covariance);
 vec Mu=Rcpp::as<vec>(obj["Mu"]), xmean=Rcpp::as<vec>(initial_mean);
 mat GM=Rcpp::as<mat>(obj["Gamma_macro"]), GS=Rcpp::as<mat>(obj["Gamma_s"]);
 vec GS0=Rcpp::as<vec>(obj["Gamma_s0"]);
 Rcpp::List pars=obj["pars"], co=obj["coefs_q"], loadings=obj["loadings"];
 vec shadow_loading=Rcpp::as<vec>(loadings["shadow_rate"]);
 Rcpp::NumericVector Af=co["A_X_for"], Ae=co["A_X_exp"], Api=co["A_X_for_pi"], sd_n=obj["s_n"];
 Rcpp::NumericMatrix Bf=co["B_X_for"], Bpi=co["B_X_pi"], SX=obj["Sigma2_X"],
   Bcum=co["B_X_cum"], Bcpi=co["B_X_cum_pi"];
 vec sigma_r=Rcpp::as<vec>(obj["sigma_r_t"]), real_offset=Rcpp::as<vec>(obj["real_yield_offset_t"]);
 const int T=Y.n_rows,K=Phi.n_rows,Ns=GS.n_rows,Nt=horizons.size(),Nn=mats_n.size(),Nr=mats_r.size(),N=Y.n_cols;
 if(N!=2+Ns+Nt+Nn+Nr || K!=6)Rcpp::stop("Unsupported observation layout in private accelerator");
 double lb=Rcpp::as<double>(pars["r_lb"]);
 vec baseR(N);
 baseR[0]=std::pow(Rcpp::as<double>(pars["sigma_inf"]),2);
 baseR[1]=std::pow(Rcpp::as<double>(pars["sigma_ptr"]),2);
 for(int j=0;j<Ns;++j)baseR[2+j]=std::pow(Rcpp::as<double>(pars[j==1?"sigma_inflation_survey_long":"sigma_inflation_survey"]),2);
 for(int j=0;j<Nt;++j)baseR[2+Ns+j]=std::pow(Rcpp::as<double>(pars[j==1?"sigma_tbill_survey_long":"sigma_tbill_survey"]),2);
 for(int j=0;j<Nn;++j)baseR[2+Ns+Nt+j]=std::pow(Rcpp::as<double>(pars["sigma_nominal_yield"]),2);
 Rcpp::List terms=obj["tbill_terms"];
 std::vector<mat> tB;std::vector<vec> tA,tSD;std::vector<int> tH;
 for(int j=0;j<terms.size();++j){Rcpp::List term=terms[j];tB.push_back(Rcpp::as<mat>(term["loading"]));
   tA.push_back(Rcpp::as<vec>(term["intercept"]));tSD.push_back(Rcpp::as<vec>(term["sd"]));tH.push_back(Rcpp::as<int>(term["horizon"]));}
 Rcpp::LogicalVector desired;if(branches.isNotNull())desired=Rcpp::LogicalVector(branches);
 int branch_index=0;std::vector<int> actual_trace;std::vector<double> shadow_trace;
 mat xp(T,K),xu(T,K);vec ll(T,arma::fill::zeros);
 Rcpp::NumericVector prior_array(K*K*T),posterior_array(K*K*T);
 prior_array.attr("dim")=Rcpp::IntegerVector::create(K,K,T);posterior_array.attr("dim")=Rcpp::IntegerVector::create(K,K,T);
 Rcpp::LogicalVector iterated(T,false);
 auto measurement=[&](const vec& state,int tt,const mat& prior)->Measurement {
   Rcpp::NumericMatrix xx(K,1);for(int k=0;k<K;++k)xx(k,0)=state[k];
   Rcpp::List yf=y_fitting_jfec_cpp(xx,Af,Bf,Ae,lb,sd_n,Api,Bpi,SX,Bcum,Bcpi);
   Rcpp::NumericMatrix yn=yf["yfit_all_n"],yr=yf["yfit_all_r"],Jn=yf["JJ_n"],Jr=yf["JJ_r"];
   double shadow=arma::dot(shadow_loading,state);int actual=shadow>lb;
   shadow_trace.push_back(1200*shadow);
   if(record)actual_trace.push_back(actual);
   int selected=actual;
   if(branches.isNotNull()){
     if(branch_index>=desired.size())Rcpp::stop("Short branch trace");selected=desired[branch_index];}
   if(record||branches.isNotNull())++branch_index;
   double correction=lb+selected*(shadow-lb)-std::max(shadow,lb);
   mat G(N,K,arma::fill::zeros);vec predicted(N);predicted.fill(arma::datum::nan);
   G.rows(0,1)=GM;predicted.head(2)=GM*state;
   G.rows(2,2+Ns-1)=GS;predicted.subvec(2,2+Ns-1)=GS0+GS*state;
   bool observed_tbill=false;
   for(int j=0;j<Nt;++j)if(std::isfinite(Y(tt,2+Ns+j)))observed_tbill=true;
   if(observed_tbill){
     vec val(terms.size());mat jac(terms.size(),K);
     for(int h=0;h<terms.size();++h){
       vec mu=tA[h]+tB[h]*state;vec prob(3);double sum=0;
       for(int r=0;r<3;++r){double z=mu[r]/tSD[h][r];prob[r]=filter_norm_cdf(z);sum+=mu[r]*prob[r]+tSD[h][r]*filter_norm_pdf(z);}
       val[h]=lb+sum/3; jac.row(h)=(prob.t()*tB[h])/3;
     }
     for(int j=0;j<Nt;++j)if(std::isfinite(Y(tt,2+Ns+j))){
       int count=0;double sum=0;arma::rowvec g(K,arma::fill::zeros);
       for(int h=0;h<terms.size();++h)if(tH[h]<horizons[j]){++count;sum+=val[h];g+=jac.row(h);}
       predicted[2+Ns+j]=sum/count;G.row(2+Ns+j)=g/count;
     }
   }
   int row=2+Ns+Nt;
   for(int j=0;j<Nn;++j){int m=mats_n[j];predicted[row+j]=yn(0,m-1)+correction/m;
     for(int k=0;k<K;++k)G(row+j,k)=Jn(m-1,k)+(selected-actual)*Bf(k,0)/m;}
   row+=Nn;
   for(int j=0;j<Nr;++j){int m=mats_r[j];predicted[row+j]=yr(0,m-1)+correction/m+real_offset[tt];
     for(int k=0;k<K;++k)G(row+j,k)=Jr(m-1,k)+(selected-actual)*Bf(k,0)/m;}
   vec R=baseR;for(int j=0;j<Nr;++j)R[row+j]=sigma_r[tt]*sigma_r[tt];
   std::vector<arma::uword> indices;
   for(int i=0;i<N;++i)if(std::isfinite(Y(tt,i))){
     double magnitude=arma::accu(arma::abs(G.row(i)));
     if(!std::isfinite(magnitude)||!std::isfinite(predicted[i]))Rcpp::stop("Non-finite measurement");
     if(magnitude>0)indices.push_back(i);
   }
   Measurement answer;arma::uvec idx(indices);answer.G=G.rows(idx);
   answer.innovation=Y.row(tt).t().eval().elem(idx)-predicted.elem(idx);
   answer.cross=prior*answer.G.t();answer.delta=answer.G*prior*answer.G.t();
   answer.delta.diag()+=R.elem(idx);answer.delta=(answer.delta+answer.delta.t())/2;
   return answer;
 };
 mat posterior=P;vec updated=xmean;
 for(int tt=0;tt<T;++tt){
   if(tt%100==0)Rcpp::checkUserInterrupt();
   if(tt>0){xmean=Mu+Phi*updated;P=Phi*posterior*Phi.t()+Sigma;P=(P+P.t())/2;}
   xp.row(tt)=xmean.t();for(int j=0;j<K;++j)for(int i=0;i<K;++i)prior_array[i+K*j+K*K*tt]=P(i,j);
   Measurement sys=measurement(xmean,tt,P);
   if(sys.innovation.n_elem==0){updated=xmean;posterior=P;}
   else {
     mat U,inv;if(!arma::chol(U,sys.delta)||!arma::inv_sympd(inv,sys.delta))Rcpp::stop("Non-SPD innovation covariance");
     ll[tt]=-.5*sys.innovation.n_elem*std::log(2*arma::datum::pi)-arma::accu(arma::log(U.diag()))-.5*arma::as_scalar(sys.innovation.t()*inv*sys.innovation);
     mat gain=sys.cross*inv;updated=xmean+gain*sys.innovation;
     double variance=arma::as_scalar(shadow_loading.t()*P*shadow_loading);
     double probability=filter_norm_cdf((lb-arma::dot(shadow_loading,xmean))/std::sqrt(variance));
     if(iterated_filter && probability>=trigger){
       iterated[tt]=true;vec expansion=updated;
       for(int step=2;step<=steps;++step){
         sys=measurement(expansion,tt,P);vec adjusted=sys.innovation+sys.G*(expansion-xmean);
         if(!arma::inv_sympd(inv,sys.delta))Rcpp::stop("Non-SPD iterated covariance");
         gain=sys.cross*inv;updated=xmean+gain*adjusted;expansion=updated;
       }
     }
     posterior=P-gain*sys.cross.t();posterior=(posterior+posterior.t())/2;
   }
   xu.row(tt)=updated.t();for(int j=0;j<K;++j)for(int i=0;i<K;++i)posterior_array[i+K*j+K*K*tt]=posterior(i,j);
 }
 if(branches.isNotNull()&&branch_index!=desired.size())Rcpp::stop("Unconsumed branch trace");
 Rcpp::LogicalVector trace(actual_trace.size());for(int i=0;i<actual_trace.size();++i)trace[i]=actual_trace[i];
 return Rcpp::List::create(Rcpp::Named("fval")=arma::accu(ll),Rcpp::Named("fval_t")=ll,
   Rcpp::Named("x_pred")=xp,Rcpp::Named("x_upd")=xu,Rcpp::Named("P_pred")=prior_array,
   Rcpp::Named("P_upd")=posterior_array,Rcpp::Named("iterated")=iterated,
   Rcpp::Named("nominal_branches")=trace,Rcpp::Named("shadow_trace_pp")=shadow_trace);
}
