// Two-factor nominal/real shadow-rate Monte Carlo. P=Q; known means.
// Rates/states are annualized percentage points; pricing uses monthly units.
#include <Rcpp.h>
using namespace Rcpp;
struct Pricing {
  double b[60][2],bp[60],a[60],ap[60],sd[60],cross[60];
  double sig[2],rho[2],means[2];
  Pricing(NumericVector th,NumericVector mu) {
    for(int k=0;k<2;k++){
      rho[k]=1/(1+exp(-th[k]));
      sig[k]=exp(th[2+k])/1200;
      means[k]=mu[k]/1200;
    }
    double c[2]= {0,0},cp=1,v=0;
    for(int j=0;j<60;j++){
      b[j][0]=pow(rho[0],j);
      b[j][1]=pow(rho[1],j);
      bp[j]=pow(rho[1],j+1);
      if(j>0){
        for(int k=0;k<2;k++)c[k]+=b[j-1][k];
        cp+=pow(rho[1],j);
      }
      a[j]=-.5*(c[0]*c[0]*sig[0]*sig[0]+c[1]*c[1]*sig[1]*sig[1]);
      ap[j]=.5*cp*cp*sig[1]*sig[1];
      cross[j]=j==0?0:c[1]*cp*sig[1]*sig[1];
      sd[j]=sqrt(v);
      v+=b[j][0]*b[j][0]*sig[0]*sig[0]+b[j][1]*b[j][1]*sig[1]*sig[1];
    }
  }
  void measure(const double*x,double*h,double g[5][2],int forced_branch=-1) {
    double cn=0,cr=0,gn[2]= {0,0},gr[2]= {0,0};
    int out=0;
    for(int j=0;j<60;j++){
      double m=means[0]+means[1]+(b[j][0]*x[0]+b[j][1]*x[1])/1200;
      double z=j?m/sd[j]:0,prob=j?R::pnorm(z,0,1,1,0):(forced_branch>=0?forced_branch:(m>0?1:0));
      double nm=m+a[j],pz=j?R::pnorm(nm/sd[j],0,1,1,0):(forced_branch>=0?forced_branch:(m>0?1:0));
      double f=j?nm*pz+sd[j]*R::dnorm(nm/sd[j],0,1,0):(forced_branch>=0?forced_branch*m:std::max(m,0.0));
      cn+=f;
      cr+=f-means[1]-bp[j]*x[1]/1200-ap[j]+prob*cross[j];
      for(int k=0;k<2;k++){
        double dg=pz*b[j][k];
        gn[k]+=dg;
        gr[k]+=dg-(k==1?bp[j]:0)+(j?R::dnorm(z,0,1,0)/sd[j]*b[j][k]*cross[j]:0);
      }
      if(j==11||j==59){
        h[out]=1200*cn/(j+1);
        h[out+1]=1200*cr/(j+1);
        for(int k=0;k<2;k++){
          g[out][k]=gn[k]/(j+1);
          g[out+1][k]=gr[k]/(j+1);
        }
        out+=2;
      }
    }
    h[4]=1200*means[1]+x[1];
    g[4][0]=0;
    g[4][1]=1;
  }
};
// [[Rcpp::export]]
NumericMatrix small_measure(NumericMatrix X,NumericVector theta,NumericVector means){
  Pricing p(theta,means);
  NumericMatrix Y(X.nrow(),5);
  for(int t=0;t<X.nrow();t++){
    double x[2]= {X(t,0),X(t,1)},h[5],g[5][2];
    p.measure(x,h,g);
    for(int j=0;j<5;j++)Y(t,j)=h[j];
  }
  return Y;
}
// [[Rcpp::export]]
NumericVector small_filter(NumericMatrix Y,NumericVector theta,NumericVector means,int iterations=2,Nullable<IntegerVector> fixed_branches=R_NilValue,Nullable<NumericVector> noise_multipliers=R_NilValue){
  Pricing p(theta,means);
  double q[2]= {pow(exp(theta[2]),2),pow(exp(theta[3]),2)};
  double prior[2]= {0,0},P[2][2]={
    {
      q[0]/(1-p.rho[0]*p.rho[0]),0
    },{
      0,q[1]/(1-p.rho[1]*p.rho[1])
    }
  };
  double noise=exp(2*theta[4]);
  NumericVector ll(Y.nrow());
  NumericVector ns(5,1.0);
  if(noise_multipliers.isNotNull()){
    ns=as<NumericVector>(noise_multipliers.get());
    if(ns.size()!=5)stop("Need five noise multipliers");
    for(int j=0;j<5;j++)if(!(ns[j]>0))stop("Positive noise multipliers required");
  }
  IntegerVector trace(Y.nrow()*iterations);
  IntegerVector active;
  bool use_fixed=fixed_branches.isNotNull();
  if(use_fixed){
    active=as<IntegerVector>(fixed_branches.get());
    if(active.size()!=trace.size())stop("Invalid fixed-branch trace");
  }
  for(int t=0;t<Y.nrow();t++){
    double expansion[2]= {prior[0],prior[1]},updated[2],U[2][2];
    for(int it=0;it<iterations;it++){
      double h[5],g[5][2];
      int branch=use_fixed?active[t*iterations+it]:-1;
      trace[t*iterations+it]=(means[0]+means[1]+expansion[0]+expansion[1]>0)?1:0;
      p.measure(expansion,h,g,branch);
      for(int k=0;k<2;k++){
        updated[k]=prior[k];
        for(int l=0;l<2;l++)U[k][l]=P[k][l];
      }
      for(int j=0;j<5;j++){
        if(!R_finite(Y(t,j)))continue;
        double r=Y(t,j)-h[j];
        for(int k=0;k<2;k++)r+=g[j][k]*(expansion[k]-updated[k]);
        double c[2]= {U[0][0]*g[j][0]+U[0][1]*g[j][1],U[1][0]*g[j][0]+U[1][1]*g[j][1]};
        double v=noise*ns[j]*ns[j]+g[j][0]*c[0]+g[j][1]*c[1];
        if(!(v>0))return NumericVector(Y.nrow(),NA_REAL);
        if(it==0)ll[t]+=-.5*(log(2*M_PI*v)+r*r/v);
        for(int k=0;k<2;k++){
          updated[k]+=c[k]/v*r;
          for(int l=0;l<2;l++)U[k][l]-=c[k]*c[l]/v;
        }
      }
      for(int k=0;k<2;k++)expansion[k]=updated[k];
    }
    for(int k=0;k<2;k++){
      prior[k]=p.rho[k]*updated[k];
      for(int l=0;l<2;l++)P[k][l]=p.rho[k]*p.rho[l]*U[k][l]+(k==l?q[k]:0);
    }
  }
  ll.attr("branches")=trace;
  return ll;
}
// Scaled unscented filter: alpha=1, beta=2, kappa=0; five sigma points.
// Optional linear map is used only for an independent exact-KF validation.
// [[Rcpp::export]]
NumericVector small_ukf(NumericMatrix Y,NumericVector theta,NumericVector means,
bool linear_check=false,Nullable<NumericVector> noise_multipliers=R_NilValue) {
  Pricing p(theta,means);
  const int N=5;
  double q[2]= {pow(exp(theta[2]),2),pow(exp(theta[3]),2)};
  double x[2]= {0,0},P[2][2]={
    {
      q[0]/(1-p.rho[0]*p.rho[0]),0
    },{
      0,q[1]/(1-p.rho[1]*p.rho[1])
    }
  };
  double noise=exp(2*theta[4]);
  NumericVector ll(Y.nrow());
  NumericVector ns(5,1.0);
  if(noise_multipliers.isNotNull()){
    ns=as<NumericVector>(noise_multipliers.get());
    if(ns.size()!=5)stop("Need five noise multipliers");
    for(int j=0;j<5;j++)if(!(ns[j]>0))stop("Positive noise multipliers required");
  }
  const double wm[5]= {0,.25,.25,.25,.25},wc[5]= {2,.25,.25,.25,.25};
  const double linear_H[5][2]={
    {
      1,1
    },{
      1,0
    },{
      .5,.5
    },{
      1,-.5
    },{
      0,1
    }
  };
  for(int t=0;t<Y.nrow();t++) {
    if(!(P[0][0]>0))return NumericVector(Y.nrow(),NA_REAL);
    double l00=sqrt(P[0][0]),l10=P[1][0]/l00,remainder=P[1][1]-l10*l10;
    if(!(remainder>0))return NumericVector(Y.nrow(),NA_REAL);
    double l11=sqrt(remainder),sp[5][2],hy[5][5],ym[5]= {0,0,0,0,0};
    for(int a=0;a<5;a++){
      sp[a][0]=x[0];
      sp[a][1]=x[1];
    }
    sp[1][0]+=sqrt(2.0)*l00;
    sp[1][1]+=sqrt(2.0)*l10;
    sp[2][1]+=sqrt(2.0)*l11;
    sp[3][0]-=sqrt(2.0)*l00;
    sp[3][1]-=sqrt(2.0)*l10;
    sp[4][1]-=sqrt(2.0)*l11;
    for(int a=0;a<5;a++) {
      double g[5][2];
      p.measure(sp[a],hy[a],g);
      if(linear_check)for(int j=0;j<N;j++)hy[a][j]=linear_H[j][0]*sp[a][0]+linear_H[j][1]*sp[a][1];
      for(int j=0;j<N;j++)ym[j]+=wm[a]*hy[a][j];
    }
    double S[5][5]={
      {
        0
      }
    },C[2][5]={
      {
        0
      }
    },L[5][5]={
      {
        0
      }
    },r[5],z[5];
    for(int j=0;j<N;j++){
      r[j]=Y(t,j)-ym[j];
      if(!R_finite(r[j]))return NumericVector(Y.nrow(),NA_REAL);
      for(int k=0;k<N;k++){
        S[j][k]=(j==k?noise*ns[j]*ns[j]:0);
        for(int a=0;a<5;a++)S[j][k]+=wc[a]*(hy[a][j]-ym[j])*(hy[a][k]-ym[k]);
      }
      for(int k=0;k<2;k++)for(int a=0;a<5;a++)C[k][j]+=wc[a]*(sp[a][k]-x[k])*(hy[a][j]-ym[j]);
    }
    double logdet=0,quad=0;
    for(int i=0;i<N;i++){
      for(int j=0;j<=i;j++){
        double v=S[i][j];
        for(int k=0;k<j;k++)v-=L[i][k]*L[j][k];
        if(i==j){
          if(!(v>0))return NumericVector(Y.nrow(),NA_REAL);
          L[i][j]=sqrt(v);
          logdet+=log(L[i][j]);
        }
        else L[i][j]=v/L[j][j];
      }
      double v=r[i];
      for(int j=0;j<i;j++)v-=L[i][j]*z[j];
      z[i]=v/L[i][i];
      quad+=z[i]*z[i];
    }
    ll[t]=-.5*N*log(2*M_PI)-logdet-.5*quad;
    double whitened_C[2][5],updated[2],U[2][2];
    for(int k=0;k<2;k++){
      updated[k]=x[k];
      for(int i=0;i<N;i++){
        double v=C[k][i];
        for(int j=0;j<i;j++)v-=L[i][j]*whitened_C[k][j];
        whitened_C[k][i]=v/L[i][i];
        updated[k]+=whitened_C[k][i]*z[i];
      }
    }
    for(int k=0;k<2;k++)for(int j=0;j<2;j++){
      U[k][j]=P[k][j];
      for(int i=0;i<N;i++)U[k][j]-=whitened_C[k][i]*whitened_C[j][i];
    }
    for(int k=0;k<2;k++){
      x[k]=p.rho[k]*updated[k];
      for(int j=0;j<2;j++)P[k][j]=p.rho[k]*p.rho[j]*U[k][j]+(k==j?q[k]:0);
    }
  }
  return ll;
}
