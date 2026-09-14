// ---- Enable C++11: ----
// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/SpecialFunctions>

//[[Rcpp::depends(RcppEigen)]] 


using namespace Rcpp ;
using namespace Eigen ;
using Eigen::RealSchur;
using Eigen::MatrixXd;

// [[Rcpp::export]]
Rcpp::List compute_schur(MatrixXd x) {
  RealSchur<MatrixXd> schur(x);
  return Rcpp::List::create(_["U"] = schur.matrixU(), 
                            _["T"] = schur.matrixT());
}

Eigen::MatrixXd Vec(const Eigen::MatrixXd & M){
  int rM = M.rows() ;
  int cM = M.cols() ;
  Eigen::MatrixXd vecM  = Eigen::MatrixXd::Zero(rM*cM, 1);
  for(int j = 1; j <= cM; j++){
    vecM.block((j-1)*rM,0,rM,1) = M.col(j-1) ;
  }
  
  return vecM ;
}

Eigen::MatrixXd Vec_1(const Eigen::MatrixXd & vecM){
  int n2 = vecM.rows() ;
  int n = sqrt(n2) ;
  
  Eigen::MatrixXd M  = Eigen::MatrixXd::Zero(n, n);
  for(int j = 1; j <= n; j++){
    M.block(0,j-1,n,1) = vecM.block((j-1)*n,0,n,1) ;
  }
  return M ;
}

Eigen::MatrixXd kronecker_cpp(const Eigen::MatrixXd & A,
                              const Eigen::MatrixXd & B){
  int rA = A.rows() ;
  int cA = A.cols() ;
  int rB = B.rows() ;
  int cB = B.cols() ;
  
  Eigen::MatrixXd M   = Eigen::MatrixXd::Zero(rA*rB, cA*cB) ;
  Eigen::MatrixXd AUX = Eigen::MatrixXd::Zero(rB, cB) ;
  
  for(int i = 1; i <= rA; i++){
    for(int j = 1; j <= cA; j++){
      AUX = Eigen::MatrixXd::Constant(rB, cB, A(i-1,j-1)) ;
      AUX = AUX.array() * B.array() ;
      M.block((i-1)*rB,(j-1)*cB,rB,cB) = AUX ;
    }
  }
  return M ;
}

Eigen::MatrixXd extract_diagonal(const Eigen::MatrixXd & M) {
  
  int n = M.rows() ;
  
  Eigen::MatrixXd diag = Eigen::MatrixXd::Zero(n,1) ;
  
  for (int i = 0; i < n; ++i) {
    diag(i,0) = M(i,i) ;
  }
  
  return diag;
}

// [[Rcpp::export]]
Eigen::MatrixXd compute_uncond_Var(const Eigen::MatrixXd & Phi,
                                   const Eigen::MatrixXd & Sigma){
  
  int n = Phi.rows() ;
  
  Eigen::MatrixXd SS = Eigen::MatrixXd::Zero(n,n) ;
  SS = Sigma * Sigma.transpose() ;
  
  Eigen::MatrixXd vecSS = Eigen::MatrixXd::Zero(n*n,1) ;
  vecSS = Vec(SS) ;
  
  Eigen::MatrixXd Idnn = Eigen::MatrixXd::Identity(n*n,n*n) ;
  
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n*n,n*n) ;
  
  M = Idnn - kronecker_cpp(Phi,Phi) ;
  M = M.inverse() ;
  
  Eigen::MatrixXd vecVAR = Eigen::MatrixXd::Zero(n*n,1) ;
  
  vecVAR = M * vecSS ;
  
  Eigen::MatrixXd VAR = Eigen::MatrixXd::Zero(n,n) ;
  VAR = Vec_1(vecVAR) ;
  
  return VAR;
}


// [[Rcpp::export]]
Rcpp::List compute_share_kth_shock(const Eigen::MatrixXd & A,
                                   const Eigen::MatrixXd & Phi_tilde,
                                   const Eigen::MatrixXd & Sigma_tilde,
                                   const int k,
                                   const int horizon){
  // This function compute the share of variance accounted for by the k^th shock.
  
  int n = Phi_tilde.rows() ;
  int nb_shocks = Sigma_tilde.cols() ;
  
  Eigen::MatrixXd Phi   = Eigen::MatrixXd::Zero(n,n) ;
  Eigen::MatrixXd A_1   = Eigen::MatrixXd::Zero(n,n) ;
  Eigen::MatrixXd Sigma = Eigen::MatrixXd::Zero(n,nb_shocks) ;
  
  A_1   = A.inverse() ;
  Phi   = A_1 * Phi_tilde ;
  Sigma = A_1 * Sigma_tilde ;
  
  Eigen::MatrixXd SS = Eigen::MatrixXd::Zero(n,n) ;
  SS = Sigma * Sigma.transpose() ;
  
  Eigen::MatrixXd Sigma_k = Eigen::MatrixXd::Zero(n,nb_shocks) ;
  Sigma_k = Sigma ;
  Eigen::MatrixXd vecnZeros = Eigen::MatrixXd::Zero(n,1) ;
  Sigma_k.block(0,k-1,n,1) = vecnZeros ;
  
  Eigen::MatrixXd SS_k = Eigen::MatrixXd::Zero(n,n) ;
  SS_k = Sigma_k * Sigma_k.transpose() ;
  
  Eigen::MatrixXd Idnn = Eigen::MatrixXd::Identity(n*n,n*n) ;
  
  Eigen::MatrixXd VAR   = Eigen::MatrixXd::Zero(n,n) ;
  Eigen::MatrixXd VAR_k = Eigen::MatrixXd::Zero(n,n) ;
  
  Eigen::MatrixXd diagVAR     = Eigen::MatrixXd::Zero(n,1) ;
  Eigen::MatrixXd diagVAR_k   = Eigen::MatrixXd::Zero(n,1) ;
  
  if(horizon>100){
    Eigen::MatrixXd vecSS   = Eigen::MatrixXd::Zero(n*n,1) ;
    Eigen::MatrixXd vecSS_k = Eigen::MatrixXd::Zero(n*n,1) ;
    vecSS   = Vec(SS) ;
    vecSS_k = Vec(SS_k) ;
    
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n*n,n*n) ;
    
    M = Idnn - kronecker_cpp(Phi,Phi) ;
    M = M.inverse() ;
    
    Eigen::MatrixXd vecVAR   = Eigen::MatrixXd::Zero(n*n,1) ;
    Eigen::MatrixXd vecVAR_k = Eigen::MatrixXd::Zero(n*n,1) ;
    
    vecVAR   = M * vecSS ;
    vecVAR_k = M * vecSS_k ;
    
    VAR   = Vec_1(vecVAR) ;
    VAR_k = Vec_1(vecVAR_k) ;
  }else{
    VAR   = SS ;
    VAR_k = SS_k ;
    for(int h = 2; h <= horizon; h++){
      VAR   = SS   + Phi * VAR   * Phi.transpose() ;
      VAR_k = SS_k + Phi * VAR_k * Phi.transpose() ;
    }
  }
  
  diagVAR   = extract_diagonal(VAR) ;
  diagVAR_k = extract_diagonal(VAR_k) ;
  
  Eigen::MatrixXd share = Eigen::MatrixXd::Zero(n,1) ;
  
  share = (diagVAR - diagVAR_k).array() * (diagVAR.cwiseInverse()).array() ;
  
  return List::create(Named("share")     = share) ;
}


Eigen::MatrixXd na_matrix(int r,
                          int c){
  // create matrix with NaNs
  Eigen::MatrixXd m(r,c) ;
  m = MatrixXd::Zero(r,c) ;
  //m(0,0) = R_NaN ;
  m.fill(R_NaN) ;
  return m ;
}

Eigen::VectorXd getEigenValues(const Eigen::MatrixXd & M) {
  SelfAdjointEigenSolver<MatrixXd> es(M);
  return es.eigenvalues();
}


Eigen::MatrixXd positiveIndicator(const Eigen::MatrixXd & M) {
  int r = M.rows() ;
  int c = M.cols() ;
  Eigen::MatrixXd Mpos = Eigen::MatrixXd::Zero(r,c) ;
  for(int i = 0; i < r; i++){
    for(int j = 1; j < c; j++){
      if(M(i,j)>=0){
        Mpos(i,j) = 1 ;
      }
    }
  }
  return Mpos;
}

Eigen::MatrixXd which_min_row(const Eigen::MatrixXd & M){
  int p = M.rows() ;
  int n = M.cols() ;
  Eigen::MatrixXd VecMins = Eigen::MatrixXd::Zero(p,1) ;
  for(int i = 0; i < p; i++){
    int aux = 0 ;
    double mini = M(i,0) ;
    double v = 0 ;
    for(int j = 1; j < n; j++){
      v = M(i,j) ;
      if( v < mini ){
        aux = j ;
        mini = v ;
      }
    }
    VecMins(i) = aux + 1;
  }
  return VecMins;
}

Eigen::MatrixXd find_closest_in_vec(const Eigen::MatrixXd & M,
                                    const Eigen::MatrixXd & vec){
  int p = M.rows() ;
  int n = M.cols() ;
  int vec_values = vec.rows() ;
  Eigen::MatrixXd MatIndicators = Eigen::MatrixXd::Zero(p,n) ;
  double smallest_dist ;
  double dist ;
  double coeff ;
  double coeff_vec ;
  for(int i = 0; i < p; i++){
    for(int j = 0; j < n; j++){
      coeff = M(i,j) ;
      coeff_vec = vec(0,0) ;
      smallest_dist = (coeff - coeff_vec) * (coeff - coeff_vec) ;
      MatIndicators(i,j) = 1 ; // may stay there if no one smaller
      for(int k = 1; k < vec_values; k++){
        coeff_vec = vec(k,0) ;
        dist = (coeff - coeff_vec) * (coeff - coeff_vec) ;
        if(dist < smallest_dist){
          MatIndicators(i,j) = k + 1;
          smallest_dist = dist ;
        }
      }
    }
  }
  return MatIndicators ;
}

Eigen::MatrixXd fill_from_indic(const Eigen::MatrixXd & MatIndic,
                                const Eigen::MatrixXd & vec){
  int p = MatIndic.rows() ;
  int n = MatIndic.cols() ;
  Eigen::MatrixXd MatOutput = Eigen::MatrixXd::Zero(p,n) ;
  for(int i = 0; i < p; i++){
    for(int j = 0; j < n; j++){
      int k = MatIndic(i,j) - 1;
      MatOutput(i,j) = vec(k,0) ;
    }
  }
  return MatOutput;
}

// [[Rcpp::export]]
Eigen::MatrixXd mean_per_row(const Eigen::MatrixXd & M){
  int p = M.rows() ;
  int n = M.cols() ;
  Eigen::MatrixXd MatOutput = Eigen::MatrixXd::Zero(p,1) ;
  for(int i = 0; i < p; i++){
    for(int j = 0; j < n; j++){
      MatOutput(i,0) = MatOutput(i,0) + 1/(n*1.0) * M(i,j) ;
    }
  }
  return MatOutput;
}

Eigen::MatrixXd pmax_cpp(const Eigen::MatrixXd & M,
                         double x){
  int p = M.rows() ;
  int n = M.cols() ;
  Eigen::MatrixXd MatOutput = M ;
  for(int i = 0; i < p; i++){
    for(int j = 0; j < n; j++){
      if(M(i,j)<x){
        MatOutput(i,j) = x ;
      }
    }
  }
  return MatOutput;
}


double max_cpp(double x,
               double y){
  double Output = x ;
  if(x<y){
    Output = y ;
  }
  return Output;
}

Eigen::VectorXd seqEigen(int start,
                         int end){
  int length = end - start + 1 ;
  Eigen::VectorXd V = Eigen::VectorXd::Zero(length) ;
  for(int i = 0; i < length; i++){
    V(i) = i + start ;
  }
  return V ;
}

// [[Rcpp::export]]
Eigen::MatrixXd Mpower(const Eigen::MatrixXd & A,
                       int   n){
  // compute A^(2^n)
  int rA = A.rows() ;
  int cA = A.cols() ;
  
  Eigen::MatrixXd M   = Eigen::MatrixXd::Zero(rA, cA) ;
  M = A ;
  
  for(int i = 1; i <= n; i++){
    M = M * M ;
  }
  return M ;
}



Eigen::MatrixXd make_indices_cpp(const Eigen::MatrixXd & M){
  int n_m = M.rows() ;
  Eigen::MatrixXd mat_indices = na_matrix((n_m-1), (n_m-1)) ;
  for(int i = 1; i < n_m; i++){
    Eigen::VectorXd V = seqEigen(2,n_m-i+1) +
      Eigen::VectorXd::Constant(n_m-i,(i-1)*n_m) ;
    mat_indices.block(i-1,i-1,n_m-i,1) = V ;
  }
  return mat_indices ;
}

Eigen::MatrixXd make_indices_V_cpp(const Eigen::MatrixXd & M){
  int n_m = M.cols() ;
  Eigen::MatrixXd mat_indices = na_matrix(n_m,n_m) ;
  for(int i = 1; i <= n_m; i++){
    Eigen::VectorXd V = seqEigen(1,n_m-i+1) +
      Eigen::VectorXd::Constant(n_m-i+1,i*n_m) ;
    mat_indices.block(i-1,i-1,n_m-i+1,1) = V ;
  }
  return mat_indices ;
}

Eigen::MatrixXd pnorm_cpp(Eigen::MatrixXd X){
  int nb_r = X.rows() ;
  int nb_c = X.cols() ;
  Eigen::MatrixXd p    = Eigen::MatrixXd::Zero(nb_r, nb_c);
  Eigen::MatrixXd Cst1 = Eigen::MatrixXd::Constant(nb_r, nb_c, 0.5) ;
  Eigen::MatrixXd Cst2 = Eigen::MatrixXd::Constant(nb_r, nb_c, -sqrt(0.5)) ;
  p = Cst1.array() * ((Cst2.array() * X.array()).erfc()).array() ;
  return p;
}

double pnorm_scalar_cpp(double x){
  double p = 0.5 * erfc(-sqrt(0.5) * x) ;
  return p ;
}

Eigen::MatrixXd dnorm_cpp(Eigen::MatrixXd X){
  double pi = 3.141592653589793 ;
  int nb_r = X.rows() ;
  int nb_c = X.cols() ;
  Eigen::MatrixXd p  = Eigen::MatrixXd::Zero(nb_r, nb_c);
  Eigen::MatrixXd Cst1 = Eigen::MatrixXd::Constant(nb_r, nb_c, - 0.5) ;
  Eigen::MatrixXd Cst2 = Eigen::MatrixXd::Constant(nb_r, nb_c, 1/(sqrt(2) * sqrt(pi))) ;
  p = Cst2.array() * (Cst1.array() * (X.array().square()).array()).exp() ;
  return p;
}

double dnorm_scalar_cpp(double x){
  double pi = 3.141592653589793 ;
  double p = 1/(sqrt(2) * sqrt(pi)) * exp(- 0.5 * (x * x)) ;
  return p;
}

Eigen::MatrixXd g(const Eigen::MatrixXd & X){
  int r = X.rows() ;
  int c = X.cols() ;
  Eigen::MatrixXd gX = Eigen::MatrixXd::Zero(r,c) ;
  gX = X.array() * pnorm_cpp(X).array() ;
  gX = gX + dnorm_cpp(X) ;
  return gX ;
}



Eigen::MatrixXd get_entries(Eigen::MatrixXd X,
                            Eigen::MatrixXd v){
  int nb_r = v.rows() ;
  int vv = 0;
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(nb_r, 1) ;
  for (int i = 1; i <= nb_r; ++i){
    vv = v(i-1);
    result(i-1) = X(vv-1) ;
  }
  return(result) ;
}

// [[Rcpp::export]]
Eigen::MatrixXd apply_cumsum(Eigen::MatrixXd M){
  int nb_c = M.cols() ;
  Eigen::MatrixXd MatriX ;
  MatriX = M ;
  for (int i = 1; i <= nb_c - 1; ++i){
    MatriX.col(i) = MatriX.col(i-1) + M.col(i) ;
  }
  return MatriX ;
}

Eigen::MatrixXd make_Gamma_cpp(Eigen::MatrixXd Psi){
  int nb_r = Psi.rows() ;
  
  Eigen::MatrixXd Gamma = Eigen::MatrixXd::Zero(nb_r, nb_r) ;
  for (int i = 0; i <= nb_r - 1; ++i){
    double aux = 0 ;
    for (int j = 0; j <= nb_r - i - 1; ++j){
      aux = aux + Psi(i + j,j) ; // max: i + j <= nb_r - 1
      Gamma(i,j) = aux ;
    }
    for (int j = nb_r - i; j <=  nb_r - 1; ++j){
      Gamma(i,j) = Gamma(i,j - 1) ;
    }
  }
  return Gamma ;
}  



// [[Rcpp::export]]
Rcpp::List compute_MU_and_GAMMA(const Eigen::MatrixXd & a,
                                const double & b,
                                const Eigen::MatrixXd & a_dot,
                                const double & b_dot,
                                const Eigen::MatrixXd & Phi_x,
                                const Eigen::MatrixXd & Mu_x,
                                const Eigen::MatrixXd & Sigma_x,
                                const int max_h){
  
  int n_x = Phi_x.rows() ;
  
  Eigen::MatrixXd Phi_x_h = Eigen::MatrixXd::Identity(n_x, n_x) ;
  
  Eigen::MatrixXd vec_1_max_h_1  = Eigen::MatrixXd::Ones(1,max_h - 1) ;
  
  Rcpp::List MU   = Rcpp::List(max_h) ;
  Eigen::MatrixXd MU_i_0      = Eigen::MatrixXd::Identity(max_h,1) ;
  Eigen::MatrixXd MU_i_1      = Eigen::MatrixXd::Identity(max_h,n_x) ;
  Eigen::MatrixXd MU_lambda_0 = Eigen::MatrixXd::Identity(max_h,1) ;
  Eigen::MatrixXd MU_lambda_1 = Eigen::MatrixXd::Identity(max_h,n_x) ;
  
  Eigen::MatrixXd Omega = Sigma_x * Sigma_x.transpose() ;
  
  Rcpp::List GAMMA          = Rcpp::List(max_h) ; // # The i^th layer is Gamma_{i,0}. In particular, the first layer is Gamma_{1,0} = Omega
  Eigen::MatrixXd gamma_h_0 = Eigen::MatrixXd::Zero(n_x, n_x) ;
  
  Eigen::MatrixXd KSI_a = na_matrix(n_x, max_h + 1) ; // the i^th column of KSI.a will be xi_{i-1}^a = t(Phi_x^i) x a
  KSI_a.col(0) = a ;
  
  Eigen::MatrixXd KSI_a_dot = na_matrix(n_x, max_h + 1) ; // the i^th column of KSI.a will be xi_{i-1}^a = t(Phi_x^i) x a
  KSI_a_dot.col(0) = a_dot ;
  
  Eigen::MatrixXd SIGMA_lambda = na_matrix(max_h,1) ;
  
  Eigen::MatrixXd all_a_h =  Eigen::MatrixXd::Zero(n_x, max_h) ;
  Eigen::MatrixXd all_b_h =  Eigen::MatrixXd::Zero(1  , max_h) ;
  
  Eigen::MatrixXd a_h = Eigen::MatrixXd::Zero(n_x, 1) ;
  Eigen::MatrixXd b_h = Eigen::MatrixXd::Zero(1,   1) ;
  
  Eigen::MatrixXd sum_Phi_x_h_1 = Eigen::MatrixXd::Zero(n_x, n_x) ;
  
  for (int h = 1; h <= max_h; h++){
    
    b_h = b_h + 0.5 * ( a_h - a_dot ).transpose() * Omega * ( a_h - a_dot ) +
      ( a_h - a_dot ).transpose() * Mu_x ;
    
    a_h = Phi_x.transpose() * ( a_h - a_dot ) ;
    
    all_a_h.col(h-1) = a_h ;
    all_b_h.col(h-1) = b_h ;
    
    sum_Phi_x_h_1 = sum_Phi_x_h_1 + Phi_x_h ;
    Phi_x_h       = Phi_x_h * Phi_x ;
    
    MU_i_0(h-1)     = b_dot + (a_dot.transpose() * (sum_Phi_x_h_1 * Mu_x))(0,0) ;
    MU_i_1.row(h-1) = a_dot.transpose() * Phi_x_h ;
    
    MU_lambda_0(h-1)     = b + (a.transpose() * (sum_Phi_x_h_1 * Mu_x))(0,0) ;
    MU_lambda_1.row(h-1) = a.transpose() * Phi_x_h ;
    
    gamma_h_0 = Omega + Phi_x * gamma_h_0 * Phi_x.transpose() ;
    
    GAMMA(h-1) = gamma_h_0 ;
    
    SIGMA_lambda(h-1) = sqrt((a.transpose() * gamma_h_0 * a)(0,0)) ;
    
    KSI_a.col(h)     = Phi_x_h * a ;
    KSI_a_dot.col(h) = Phi_x_h * a_dot ;
    
  }
  
  Eigen::MatrixXd PSI_a       = Eigen::MatrixXd::Zero(max_h + 1, max_h + 1) ;
  Eigen::MatrixXd PSI_a_dot   = Eigen::MatrixXd::Zero(max_h + 1, max_h + 1) ;
  Eigen::MatrixXd PSI_a_dot_a = Eigen::MatrixXd::Zero(max_h + 1, max_h + 1) ;
  Eigen::MatrixXd PSI_a_a_dot = Eigen::MatrixXd::Zero(max_h + 1, max_h + 1) ;
  Eigen::MatrixXd PSI_dot_dot = Eigen::MatrixXd::Zero(max_h + 1, max_h + 1) ;
  PSI_a                       = KSI_a.transpose()     * Omega * KSI_a ;
  PSI_a_dot                   = KSI_a_dot.transpose() * Omega * KSI_a_dot ;
  PSI_a_dot_a                 = KSI_a_dot.transpose() * Omega * KSI_a ;
  PSI_a_a_dot                 = KSI_a.transpose()     * Omega * KSI_a_dot ;
  PSI_dot_dot                 = (KSI_a + KSI_a_dot).transpose() * Omega * (KSI_a + KSI_a_dot) ;
  
  Eigen::MatrixXd mat_GAMMA           = Eigen::MatrixXd::Zero(max_h + 1, max_h + 1) ;
  Eigen::MatrixXd mat_GAMMA_dot       = Eigen::MatrixXd::Zero(max_h + 1, max_h + 1) ;
  Eigen::MatrixXd mat_GAMMA_left_dot  = Eigen::MatrixXd::Zero(max_h + 1, max_h + 1) ;
  Eigen::MatrixXd mat_GAMMA_right_dot = Eigen::MatrixXd::Zero(max_h + 1, max_h + 1) ;
  Eigen::MatrixXd mat_GAMMA_dot_dot   = Eigen::MatrixXd::Zero(max_h + 1, max_h + 1) ;
  mat_GAMMA           = make_Gamma_cpp(PSI_a) ;
  mat_GAMMA_dot       = make_Gamma_cpp(PSI_a_dot) ;
  mat_GAMMA_left_dot  = make_Gamma_cpp(PSI_a_dot_a) ;
  mat_GAMMA_right_dot = make_Gamma_cpp(PSI_a_a_dot) ;
  mat_GAMMA_dot_dot   = make_Gamma_cpp(PSI_dot_dot) ;
  
  Eigen::MatrixXd  mat_indices = Eigen::MatrixXd::Zero(max_h  , max_h) ;
  mat_indices = make_indices_cpp(PSI_a_dot) ;
  
  return List::create(Named("MU_i_0")              = MU_i_0,
                      Named("MU_i_1")              = MU_i_1,
                      Named("MU_lambda_0")         = MU_lambda_0,
                      Named("MU_lambda_1")         = MU_lambda_1,
                      Named("SIGMA_lambda")        = SIGMA_lambda,
                      Named("all.a.h")             = all_a_h,
                      Named("all.b.h")             = all_b_h,
                      Named("mat_GAMMA")           = mat_GAMMA,
                      Named("mat_GAMMA_dot")       = mat_GAMMA_dot,
                      Named("mat_GAMMA_left_dot")  = mat_GAMMA_left_dot,
                      Named("mat_GAMMA_right_dot") = mat_GAMMA_right_dot,
                      Named("mat_GAMMA_dot_dot")   = mat_GAMMA_dot_dot,
                      Named("mat_indices")         = mat_indices
  ) ;
}



// [[Rcpp::export]]
Rcpp::List compute_condit_Exp_aux_cpp(const Eigen::MatrixXd & MU_i_0,
                                      const Eigen::MatrixXd & MU_i_1,
                                      const Eigen::MatrixXd & MU_lambda_0,
                                      const Eigen::MatrixXd & MU_lambda_1,
                                      const Eigen::MatrixXd & SIGMA_lambda,
                                      const Eigen::MatrixXd & mat_GAMMA,
                                      const Eigen::MatrixXd & mat_GAMMA_dot,
                                      const Eigen::MatrixXd & mat_GAMMA_left_dot,
                                      const Eigen::MatrixXd & mat_GAMMA_right_dot,
                                      const Eigen::MatrixXd & mat_GAMMA_dot_dot,
                                      const Eigen::MatrixXd & mat_indices,
                                      const Eigen::MatrixXd & X
){
  
  int n_x   = MU_i_1.cols() ;
  int max_h = MU_i_1.rows() ;
  int T     = X.rows() ;
  
  Eigen::MatrixXd vec_1_T        = Eigen::MatrixXd::Ones(T,1) ;
  Eigen::MatrixXd vec_1_max_h_1  = Eigen::MatrixXd::Ones(1,max_h - 1) ;
  
  Eigen::MatrixXd mat_SIGMA_lambda = Eigen::MatrixXd::Ones(T,max_h) ;
  mat_SIGMA_lambda = vec_1_T * SIGMA_lambda.transpose() ;
  
  Eigen::MatrixXd MU_i      = Eigen::MatrixXd::Ones(T,max_h) ;
  Eigen::MatrixXd MU_lambda = Eigen::MatrixXd::Ones(T,max_h) ;
  
  MU_i      = vec_1_T * MU_i_0.transpose()      + X * MU_i_1.transpose() ;
  MU_lambda = vec_1_T * MU_lambda_0.transpose() + X * MU_lambda_1.transpose() ;
  
  Eigen::MatrixXd  AUX_lambda =  Eigen::MatrixXd::Zero(T, max_h) ;
  AUX_lambda = MU_lambda.array() * (mat_SIGMA_lambda.cwiseInverse()).array() ;
  
  Eigen::MatrixXd P = Eigen::MatrixXd::Zero(T, max_h) ;
  P                 = pnorm_cpp(AUX_lambda) ;
  
  Eigen::MatrixXd MU_0 = Eigen::MatrixXd::Zero(n_x, T) ;
  MU_0 = X.transpose() ;
  
  Eigen::MatrixXd MU_i_lagged = Eigen::MatrixXd::Zero(T, max_h) ;
  
  MU_i_lagged = MU_i ;
  
  Eigen::MatrixXd F_n_1_n      = Eigen::MatrixXd::Zero(T, max_h) ;
  Eigen::MatrixXd F_n_1_n_star = Eigen::MatrixXd::Zero(T, max_h) ;
  
  Eigen::MatrixXd AUX_F_n_1     = Eigen::MatrixXd::Zero(T, max_h) ;
  Eigen::MatrixXd AUX_F_n_2     = Eigen::MatrixXd::Zero(T, max_h) ;
  AUX_F_n_1 = P.array() * MU_lambda.array() ;
  AUX_F_n_2 = dnorm_cpp(-AUX_lambda).array() * mat_SIGMA_lambda.array() ;
  
  F_n_1_n = MU_i_lagged + AUX_F_n_1 + AUX_F_n_2 ;
  Eigen::MatrixXd Const_05 = Eigen::MatrixXd::Constant(T,max_h-1,0.5) ;
  Eigen::MatrixXd Const_1  = Eigen::MatrixXd::Constant(T,max_h-1,1) ;
  Eigen::MatrixXd AUX      = Eigen::MatrixXd::Zero(T,max_h-1) ;
  AUX = (P.block(0,1,T,max_h-1).array() *
    (vec_1_T * mat_GAMMA_dot_dot.block(0,1,1,max_h-1)).array()) +
    ((Const_1 - P.block(0,1,T,max_h-1)).array() *
    (vec_1_T * mat_GAMMA_dot.block(0,1,1,max_h-1)).array()) ;
  AUX = Const_05.array() * AUX.array() ;
  F_n_1_n.block(0,1,T,max_h-1) = F_n_1_n.block(0,1,T,max_h-1) - AUX ;
  
  F_n_1_n_star = MU_i_lagged ;
  AUX = (P.block(0,0,T,max_h - 1).array() * MU_lambda.block(0,0,T,max_h - 1).array()) +
    (dnorm_cpp(-AUX_lambda.block(0,0,T,max_h - 1)).array() *
    mat_SIGMA_lambda.block(0,0,T,max_h - 1).array()) ;
  F_n_1_n_star.block(0,1,T,max_h-1) = F_n_1_n_star.block(0,1,T,max_h-1) + AUX ;
  
  AUX = (P.block(0,1,T,max_h-1).array() *
    (vec_1_T * mat_GAMMA_dot.block(0,1,1,max_h-1)).array()) +
    ((Const_1 - P.block(0,1,T,max_h-1)).array() *
    (vec_1_T * mat_GAMMA_dot.block(0,1,1,max_h-1)).array()) ;
  AUX = Const_05.array() * AUX.array() ;
  F_n_1_n_star.block(0,1,T,max_h-1) = F_n_1_n_star.block(0,1,T,max_h-1) - AUX ;
  
  AUX = P.block(0,1,T,max_h-1).array() *
    (vec_1_T * mat_GAMMA_right_dot.block(0,1,1,max_h-1)).array() ;
  F_n_1_n_star.block(0,1,T,max_h-1) = F_n_1_n_star.block(0,1,T,max_h-1) - AUX ;
  
  AUX = P.block(0,1,T,max_h-1).array() *
    (vec_1_T * mat_GAMMA.block(0,1,1,max_h-1)).array() ;
  AUX = Const_05.array() * AUX.array() ;
  F_n_1_n_star.block(0,1,T,max_h-1) = F_n_1_n_star.block(0,1,T,max_h-1) - AUX ;
  
  AUX = - (vec_1_T * mat_GAMMA_dot.block(max_h,1,1,max_h-1)).array() +
    ((P.block(0,0,T,1) * vec_1_max_h_1).array() *
    (vec_1_T * mat_GAMMA_right_dot.block(max_h - 1,0,1,max_h - 1)).array()) ;
  F_n_1_n_star.block(0,1,T,max_h-1) = F_n_1_n_star.block(0,1,T,max_h-1) + AUX ;
  
  Eigen::MatrixXd SUM_GAMMA_P_n_1_n      =  Eigen::MatrixXd::Zero(T  , max_h) ;
  Eigen::MatrixXd SUM_GAMMA_P_n_1_n_star =  Eigen::MatrixXd::Zero(T  , max_h) ;
  
  Eigen::MatrixXd aux_indices ;
  Eigen::MatrixXd aux_indices_plus_one ;
  Eigen::MatrixXd AUX1 ;
  Eigen::MatrixXd AUX2 ;
  Eigen::MatrixXd AUX10 ;
  Eigen::MatrixXd AUX20 ;
  Eigen::MatrixXd AUX30 ;
  Eigen::MatrixXd AUX40 ;
  Eigen::MatrixXd AUX_SUM_GAMMA_P_n_1_n ;
  Eigen::MatrixXd AUX_SUM_GAMMA_P_n_1_n_star ;
  
  for (int n = 3; n <= max_h; n++){
    
    AUX_SUM_GAMMA_P_n_1_n = Eigen::MatrixXd::Zero(T, n-2) ;
    
    aux_indices          = mat_indices.block(n-3,0,1,n-2).transpose() ;
    aux_indices_plus_one = mat_indices.block(n-3,0,1,n-2).transpose().array() + 1 ;
    
    AUX1    = get_entries(mat_GAMMA_dot_dot,aux_indices_plus_one) ;
    AUX2    = get_entries(mat_GAMMA_dot,aux_indices_plus_one) ;
    Const_1  = Eigen::MatrixXd::Constant(T,n-2,1) ;
    AUX_SUM_GAMMA_P_n_1_n =
      (P.block(0,0,T,n-2).array()) *
      (vec_1_T * AUX1.transpose()).array() +
      (Const_1 - P.block(0,0,T,n-2)).array() *
      (vec_1_T * AUX2.transpose()).array() ;
    
    AUX10    = get_entries(mat_GAMMA_dot,aux_indices_plus_one) ;
    AUX20    = get_entries(mat_GAMMA_right_dot,aux_indices) ;
    AUX30    = get_entries(mat_GAMMA_left_dot,aux_indices_plus_one) ;
    AUX40    = get_entries(mat_GAMMA,aux_indices) ;
    
    AUX_SUM_GAMMA_P_n_1_n_star =
      (P.block(0,0,T,n-2).array()) *
      (vec_1_T * (AUX10 + AUX20 + AUX30 + AUX40).transpose()).array() +
      (Const_1 - P.block(0,0,T,n-2)).array() *
      (vec_1_T * AUX10.transpose()).array() ;
    
    SUM_GAMMA_P_n_1_n.col(n-1)      = AUX_SUM_GAMMA_P_n_1_n * (Eigen::MatrixXd::Ones(n-2,1)) ;
    SUM_GAMMA_P_n_1_n_star.col(n-1) = AUX_SUM_GAMMA_P_n_1_n_star * (Eigen::MatrixXd::Ones(n-2,1)) ;
  }
  
  F_n_1_n      = F_n_1_n      - SUM_GAMMA_P_n_1_n ;
  F_n_1_n_star = F_n_1_n_star - SUM_GAMMA_P_n_1_n_star ;
  
  Eigen::MatrixXd cumsum_f_n_1_n      = Eigen::MatrixXd::Zero(T, max_h) ;
  Eigen::MatrixXd cumsum_f_n_1_n_star = Eigen::MatrixXd::Zero(T, max_h) ;
  
  cumsum_f_n_1_n      = - apply_cumsum(F_n_1_n) ;
  cumsum_f_n_1_n_star = - apply_cumsum(F_n_1_n_star) ;
  
  // Declare outputs:
  Eigen::MatrixXd E_n      = Eigen::MatrixXd::Zero(T, max_h) ;
  Eigen::MatrixXd E_n_star = Eigen::MatrixXd::Zero(T, max_h) ;
  
  E_n      = cumsum_f_n_1_n.array().exp() ;
  E_n_star = cumsum_f_n_1_n_star.array().exp() ;
  
  return List::create(Named("E.n")        = E_n,
                      Named("E.n.star")   = E_n_star) ;
}



// [[Rcpp::export]]
Rcpp::List compute_condit_Exp_cpp(const Eigen::MatrixXd & a,
                                  const double & b,
                                  const Eigen::MatrixXd & a_dot,
                                  const double & b_dot,
                                  const Eigen::MatrixXd & Phi_x,
                                  const Eigen::MatrixXd & Mu_x,
                                  const Eigen::MatrixXd & Sigma_x,
                                  const Eigen::MatrixXd & X,
                                  const int max_h){
  
  Rcpp::List RES_AUX = compute_MU_and_GAMMA(a,b,a_dot,b_dot,
                                            Phi_x,Mu_x,Sigma_x,max_h) ;
  
  Eigen::MatrixXd MU_i_0               = RES_AUX("MU_i_0") ;
  Eigen::MatrixXd MU_i_1               = RES_AUX("MU_i_1") ;
  Eigen::MatrixXd MU_lambda_0          = RES_AUX("MU_lambda_0") ;
  Eigen::MatrixXd MU_lambda_1          = RES_AUX("MU_lambda_1") ;
  Eigen::MatrixXd SIGMA_lambda         = RES_AUX("SIGMA_lambda") ;
  Eigen::MatrixXd mat_GAMMA            = RES_AUX("mat_GAMMA") ;
  Eigen::MatrixXd mat_GAMMA_dot        = RES_AUX("mat_GAMMA_dot") ;
  Eigen::MatrixXd mat_GAMMA_left_dot   = RES_AUX("mat_GAMMA_left_dot") ;
  Eigen::MatrixXd mat_GAMMA_right_dot  = RES_AUX("mat_GAMMA_right_dot") ;
  Eigen::MatrixXd mat_GAMMA_dot_dot    = RES_AUX("mat_GAMMA_dot_dot") ;
  Eigen::MatrixXd mat_indices          = RES_AUX("mat_indices") ;
  
  Rcpp::List RES_AUX_AUX = compute_condit_Exp_aux_cpp(MU_i_0,MU_i_1,
                                                      MU_lambda_0,MU_lambda_1,
                                                      SIGMA_lambda,
                                                      mat_GAMMA,mat_GAMMA_dot,
                                                      mat_GAMMA_left_dot,mat_GAMMA_right_dot,
                                                      mat_GAMMA_dot_dot,mat_indices,X) ;
  
  Eigen::MatrixXd E_n               = RES_AUX_AUX("E.n") ;
  Eigen::MatrixXd E_n_star          = RES_AUX_AUX("E.n.star") ;
  Eigen::MatrixXd all_a_h           = RES_AUX("all.a.h") ;
  Eigen::MatrixXd all_b_h           = RES_AUX("all.b.h") ;
  
  return List::create(Named("E.n")        = E_n,
                      Named("E.n.star")   = E_n_star,
                      Named("all.a.h")    = all_a_h,
                      Named("all.b.h")    = all_b_h) ;
}






// 
// // [[Rcpp::export]]
// Rcpp::List apply_Lemma_cp_OLDp(const Eigen::MatrixXd & mu,
//                                const Eigen::MatrixXd & Phi,
//                                const Eigen::MatrixXd & Sigma,
//                                const double & a,
//                                const Eigen::MatrixXd & b,
//                                const Eigen::MatrixXd & c,
//                                const Eigen::MatrixXd & X,
//                                const int max_h){
//   
//   int n_x = Phi.cols() ;
//   int T   = X.rows() ;
//   
//   Eigen::MatrixXd SS = Eigen::MatrixXd::Zero(n_x,n_x) ;
//   SS = Sigma * Sigma.transpose() ;
//   
//   Eigen::MatrixXd vecSS = Eigen::MatrixXd::Zero(n_x*n_x,1) ;
//   vecSS = Vec(SS) ;
//   
//   // Initialization ------------------------------------------------------------
//   Eigen::MatrixXd all_a_n       = Eigen::MatrixXd::Zero(max_h,1) ;
//   Eigen::MatrixXd all_a_bar_n   = Eigen::MatrixXd::Zero(max_h,1) ;
//   Eigen::MatrixXd all_c_dot_n   = Eigen::MatrixXd::Zero(max_h,1) ;
//   Eigen::MatrixXd all_sigma_n   = Eigen::MatrixXd::Zero(max_h,1) ;
//   Eigen::MatrixXd all_b_n       = Eigen::MatrixXd::Zero(max_h,n_x) ;
//   Eigen::MatrixXd all_c_n       = Eigen::MatrixXd::Zero(max_h,n_x) ;
//   Eigen::MatrixXd all_b_tilde_n = Eigen::MatrixXd::Zero(max_h,n_x) ;
//   Eigen::MatrixXd all_c_tilde_n = Eigen::MatrixXd::Zero(max_h,n_x) ;
//   
//   Eigen::MatrixXd Phi_low_n   = Eigen::MatrixXd::Zero(n_x,n_x) ;
//   Eigen::MatrixXd Phi_low_n_1 = Eigen::MatrixXd::Zero(n_x,n_x) ;
//   Eigen::MatrixXd Phi_exp_n   = Eigen::MatrixXd::Zero(n_x,n_x) ;
//   Eigen::MatrixXd Phi_exp_n_1 = Eigen::MatrixXd::Identity(n_x,n_x) ;
//   double sigma2_n_1 = 0 ;
//   double sigma2_n   = 0 ;
//   
//   double a_n ;
//   double a_bar_n ;
//   double c_dot_n ;
//   Eigen::MatrixXd b_n   = Eigen::MatrixXd::Zero(n_x,1) ;
//   Eigen::MatrixXd c_n   = Eigen::MatrixXd::Zero(n_x,1) ;
//   
//   // Loop on maturities --------------------------------------------------------
//   
//   for (int i = 1; i <= max_h; i++){
//     
//     Phi_low_n = Phi_low_n_1 + Phi_exp_n_1 ;
//     Phi_exp_n = Phi * Phi_exp_n_1 ;
//     
//     sigma2_n = sigma2_n_1 + 
//       (b.transpose() * Phi_exp_n_1 * SS * Phi_exp_n_1.transpose() * b)(0,0) ;
//     
//     
//     a_bar_n = a + (b.transpose() * Phi_low_n * mu)(0,0) ;
//     a_n     = a_bar_n - .5 * (b.transpose() * Phi_low_n * SS * Phi_low_n.transpose() * b)(0,0) ;
//     
//     c_dot_n = (c.transpose() * Phi_low_n * mu)(0,0) + 
//                 .5 * (c.transpose() * Phi_low_n * SS * Phi_low_n.transpose() * c)(0,0) ;
//     
//     
//     // Store results -----------------------------------------------------------
//     all_a_n(i-1,0)     = a_n ;
//     all_a_bar_n(i-1,0) = a_bar_n ;
//     all_b_n.row(i-1)   = b.transpose() * Phi_exp_n ;
//     all_c_dot_n(i-1,0) = c_dot_n ;
//     all_c_n.row(i-1)   = c.transpose() * Phi_exp_n ;
//     all_sigma_n(i-1,0) = sqrt(sigma2_n) ;
//     
//     all_b_tilde_n.row(i-1) = b.transpose() * Phi_low_n ;
//     all_c_tilde_n.row(i-1) = c.transpose() * Phi_low_n ;
//     
//     // For next iteration ------------------------------------------------------
//     Phi_low_n_1 = Phi_low_n ;
//     Phi_exp_n_1 = Phi_exp_n ;
//     sigma2_n_1  = sigma2_n ;
//   }
//   
//   Eigen::MatrixXd vec1N   = Eigen::MatrixXd::Ones(max_h,1) ;
//   Eigen::MatrixXd vec0N   = Eigen::MatrixXd::Zero(max_h,1) ;
//   Eigen::MatrixXd vec1T   = Eigen::MatrixXd::Ones(T,1) ;
//   Eigen::MatrixXd vec1x   = Eigen::MatrixXd::Ones(n_x,1) ;
//   Eigen::MatrixXd vec1xx  = Eigen::MatrixXd::Ones(n_x*n_x,1) ;
//   
//   Eigen::MatrixXd aux_an_bnX    = Eigen::MatrixXd::Zero(T,max_h) ;
//   Eigen::MatrixXd aux_abarn_bnX = Eigen::MatrixXd::Zero(T,max_h) ;
//   Eigen::MatrixXd aux_s         = Eigen::MatrixXd::Zero(T,max_h) ;
//   Eigen::MatrixXd aux_s_bar     = Eigen::MatrixXd::Zero(T,max_h) ;
//   Eigen::MatrixXd aux_sigma     = Eigen::MatrixXd::Zero(T,max_h) ;
//   Eigen::MatrixXd aux_pi        = Eigen::MatrixXd::Zero(T,max_h) ;
//   
//   Eigen::MatrixXd all_nu_n_c_equal_0 = Eigen::MatrixXd::Zero(T,max_h) ;
//   Eigen::MatrixXd all_nu_n_c         = Eigen::MatrixXd::Zero(T,max_h) ;
//   
//   Eigen::MatrixXd F_c_equal_0   = Eigen::MatrixXd::Zero(T,max_h) ;
//   Eigen::MatrixXd F_c           = Eigen::MatrixXd::Zero(T,max_h) ;
//   Eigen::MatrixXd Phi_aux_s_bar = Eigen::MatrixXd::Zero(T,max_h) ;
//   Eigen::MatrixXd phi_aux_s_bar = Eigen::MatrixXd::Zero(T,max_h) ;
//   Eigen::MatrixXd Phi_aux_s     = Eigen::MatrixXd::Zero(T,max_h) ;
//   
//   Eigen::MatrixXd aux_SS_aux    = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
//   Eigen::MatrixXd aux_SS        = Eigen::MatrixXd::Zero(max_h,1) ;
//   
//   aux_an_bnX    = vec1T * all_a_n.transpose() + X * all_b_n.transpose() ;
//   aux_abarn_bnX = vec1T * all_a_bar_n.transpose() + X * all_b_n.transpose() ;
//   aux_s         = aux_an_bnX.array() * ((vec1T * all_sigma_n.transpose()).cwiseInverse()).array() ;
//   aux_s_bar     = aux_abarn_bnX.array() * ((vec1T * all_sigma_n.transpose()).cwiseInverse()).array() ;
//   
//   aux_sigma = vec1T * all_sigma_n.transpose() ;
//   aux_pi    = vec1T * all_c_dot_n.transpose() + X * all_c_n.transpose() ;
//   
//   aux_SS_aux = kronecker_cpp(all_b_tilde_n,vec1x.transpose()).array() * 
//     kronecker_cpp(vec1x.transpose(),all_c_tilde_n).array() * (vec1N * vecSS.transpose()).array() ;
//   aux_SS = aux_SS_aux * vec1xx ;
//   
//   F_c_equal_0 = aux_sigma.array() * g(aux_s).array() ;
//   
//   Phi_aux_s_bar = pnorm_cpp(aux_s_bar) ;
//   Phi_aux_s     = pnorm_cpp(aux_s) ;
//   
//   Eigen::MatrixXd auxaux    = Eigen::MatrixXd::Zero(max_h,1) ;
//   auxaux = aux_SS.array()*(all_sigma_n.cwiseInverse()).array() ;
//   phi_aux_s_bar = dnorm_cpp(aux_s_bar).array() * (vec1T * auxaux.transpose()).array() ;
//   
//   Eigen::MatrixXd indic_sigma0 = Eigen::MatrixXd::Zero(max_h,1) ;
//   for (int i = 0; i < max_h; i++){
//     if(all_sigma_n(i,0)==0){
//       indic_sigma0(i) = 1 ;
//       
//       F_c_equal_0.col(i)   = pmax_cpp(aux_an_bnX.col(i), 0) ;
//       Phi_aux_s_bar.col(i) = positiveIndicator(aux_s_bar.col(i)) ;
//       Phi_aux_s.col(i)     = positiveIndicator(aux_s.col(i)) ;
//       phi_aux_s_bar.col(i) = vec0N ;
//     }
//   }
//   
//   F_c = Phi_aux_s_bar.array() * (vec1T * aux_SS.transpose()).array() ;
//   F_c = F_c + F_c_equal_0 - aux_pi ;
//   
//   all_nu_n_c_equal_0 = Phi_aux_s ;
//   all_nu_n_c         = Phi_aux_s + phi_aux_s_bar ;
//   
//   return List::create(Named("F_c_equal_0")        = F_c_equal_0,
//                       Named("F_c")                = F_c,
//                       Named("all_b_n")            = all_b_n,
//                       Named("all_c_n")            = all_c_n,
//                       Named("all_nu_n_c")         = all_nu_n_c,
//                       Named("all_nu_n_c_equal_0") = all_nu_n_c_equal_0) ;
// }


// [[Rcpp::export]]
Rcpp::List apply_Lemma_cpp_aux(const Eigen::MatrixXd & mu,
                               const Eigen::MatrixXd & Phi,
                               const Eigen::MatrixXd & Sigma,
                               const double          & a,
                               const Eigen::MatrixXd & b,
                               const Eigen::MatrixXd & c,
                               const int max_h){
  
  int n_x = Phi.cols() ;
  
  Eigen::MatrixXd SS = Eigen::MatrixXd::Zero(n_x,n_x) ;
  SS = Sigma * Sigma.transpose() ;
  
  Eigen::MatrixXd vecSS = Eigen::MatrixXd::Zero(n_x*n_x,1) ;
  vecSS = Vec(SS) ;
  
  // Initialization ------------------------------------------------------------
  Eigen::MatrixXd all_a_n       = Eigen::MatrixXd::Zero(max_h,1) ;
  Eigen::MatrixXd all_a_bar_n   = Eigen::MatrixXd::Zero(max_h,1) ;
  Eigen::MatrixXd all_c_dot_n   = Eigen::MatrixXd::Zero(max_h,1) ;
  Eigen::MatrixXd all_sigma_n   = Eigen::MatrixXd::Zero(max_h,1) ;
  Eigen::MatrixXd all_b_n       = Eigen::MatrixXd::Zero(max_h,n_x) ;
  Eigen::MatrixXd all_c_n       = Eigen::MatrixXd::Zero(max_h,n_x) ;
  Eigen::MatrixXd all_b_tilde_n = Eigen::MatrixXd::Zero(max_h,n_x) ;
  Eigen::MatrixXd all_c_tilde_n = Eigen::MatrixXd::Zero(max_h,n_x) ;
  
  Eigen::MatrixXd Phi_low_n   = Eigen::MatrixXd::Zero(n_x,n_x) ;
  Eigen::MatrixXd Phi_low_n_1 = Eigen::MatrixXd::Zero(n_x,n_x) ;
  Eigen::MatrixXd Phi_exp_n   = Eigen::MatrixXd::Zero(n_x,n_x) ;
  Eigen::MatrixXd Phi_exp_n_1 = Eigen::MatrixXd::Identity(n_x,n_x) ;
  double sigma2_n_1 = 0 ;
  double sigma2_n   = 0 ;
  
  double a_n ;
  double a_bar_n ;
  double c_dot_n ;
  Eigen::MatrixXd b_n   = Eigen::MatrixXd::Zero(n_x,1) ;
  Eigen::MatrixXd c_n   = Eigen::MatrixXd::Zero(n_x,1) ;
  
  // Loop on maturities --------------------------------------------------------
  
  for (int i = 1; i <= max_h; i++){
    
    Phi_low_n = Phi_low_n_1 + Phi_exp_n_1 ;
    Phi_exp_n = Phi * Phi_exp_n_1 ;
    
    sigma2_n = sigma2_n_1 + 
      (b.transpose() * Phi_exp_n_1 * SS * Phi_exp_n_1.transpose() * b)(0,0) ;
    
    a_bar_n = a + (b.transpose() * Phi_low_n * mu)(0,0) ;
    a_n     = a_bar_n - .5 * (b.transpose() * Phi_low_n * SS * Phi_low_n.transpose() * b)(0,0) ;
    
    c_dot_n = (c.transpose() * Phi_low_n * mu)(0,0) + 
                .5 * (c.transpose() * Phi_low_n * SS * Phi_low_n.transpose() * c)(0,0) ;
    
    // Store results -----------------------------------------------------------
    all_a_n(i-1,0)     = a_n ;
    all_a_bar_n(i-1,0) = a_bar_n ;
    all_b_n.row(i-1)   = b.transpose() * Phi_exp_n ;
    all_c_dot_n(i-1,0) = c_dot_n ;
    all_c_n.row(i-1)   = c.transpose() * Phi_exp_n ;
    all_sigma_n(i-1,0) = sqrt(sigma2_n) ;
    
    all_b_tilde_n.row(i-1) = b.transpose() * Phi_low_n ;
    all_c_tilde_n.row(i-1) = c.transpose() * Phi_low_n ;
    
    // For next iteration ------------------------------------------------------
    Phi_low_n_1 = Phi_low_n ;
    Phi_exp_n_1 = Phi_exp_n ;
    sigma2_n_1  = sigma2_n ;
  }
  
  Eigen::MatrixXd aux_SS_aux    = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
  Eigen::MatrixXd aux_SS        = Eigen::MatrixXd::Zero(max_h,1) ;
  
  Eigen::MatrixXd vec1N   = Eigen::MatrixXd::Ones(max_h,1) ;
  Eigen::MatrixXd vec0N   = Eigen::MatrixXd::Zero(max_h,1) ;
  Eigen::MatrixXd vec1x   = Eigen::MatrixXd::Ones(n_x,1) ;
  Eigen::MatrixXd vec1xx  = Eigen::MatrixXd::Ones(n_x*n_x,1) ;
  
  aux_SS_aux = kronecker_cpp(all_b_tilde_n,vec1x.transpose()).array() * 
    kronecker_cpp(vec1x.transpose(),all_c_tilde_n).array() * (vec1N * vecSS.transpose()).array() ;
  aux_SS = aux_SS_aux * vec1xx ;
  // aux_SS is the vector of b'Phi_n S S' Phi_n' c for different n in {1,...,N}.
  
  return List::create(Named("all_a_n")        = all_a_n,
                      Named("all_a_bar_n")    = all_a_bar_n,
                      Named("all_c_dot_n")    = all_c_dot_n,
                      Named("all_sigma_n")    = all_sigma_n,
                      Named("all_b_n")        = all_b_n,
                      Named("all_c_n")        = all_c_n,
                      Named("all_b_tilde_n")  = all_b_tilde_n,
                      Named("all_c_tilde_n")  = all_c_tilde_n,
                      Named("aux_SS")         = aux_SS,
                      Named("vecSS")          = vecSS) ;
}



// [[Rcpp::export]]
Rcpp::List apply_Lemma_cpp_from_RES_OLD(const Rcpp::List      & RES_aux,
                                        const Eigen::MatrixXd & X,
                                        const double i_bar,
                                        const int max_h){
  
  Eigen::MatrixXd all_a_n       = RES_aux("all_a_n") ;
  Eigen::MatrixXd all_a_bar_n   = RES_aux("all_a_bar_n") ;
  Eigen::MatrixXd all_c_dot_n   = RES_aux("all_c_dot_n") ;
  Eigen::MatrixXd all_sigma_n   = RES_aux("all_sigma_n") ;
  Eigen::MatrixXd all_b_n       = RES_aux("all_b_n") ;
  Eigen::MatrixXd all_c_n       = RES_aux("all_c_n") ;
  Eigen::MatrixXd all_b_tilde_n = RES_aux("all_b_tilde_n") ;
  Eigen::MatrixXd all_c_tilde_n = RES_aux("all_c_tilde_n") ;
  Eigen::MatrixXd aux_SS        = RES_aux("aux_SS") ;
  Eigen::MatrixXd vecSS         = RES_aux("vecSS") ;
  
  int n_x = all_b_n.cols() ;
  int T   = X.rows() ;
  
  Eigen::MatrixXd vec1N   = Eigen::MatrixXd::Ones(max_h,1) ;
  Eigen::MatrixXd vec0N   = Eigen::MatrixXd::Zero(max_h,1) ;
  Eigen::MatrixXd vec1T   = Eigen::MatrixXd::Ones(T,1) ;
  Eigen::MatrixXd vec1x   = Eigen::MatrixXd::Ones(n_x,1) ;
  Eigen::MatrixXd vec1xx  = Eigen::MatrixXd::Ones(n_x*n_x,1) ;
  
  Eigen::MatrixXd aux_an_bnX    = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd aux_abarn_bnX = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd aux_s         = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd aux_s_bar     = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd aux_sigma     = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd aux_pi        = Eigen::MatrixXd::Zero(T,max_h) ;
  
  // For derivative:
  Eigen::MatrixXd all_nu_n_c_equal_0 = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd all_nu_n_c         = Eigen::MatrixXd::Zero(T,max_h) ;
  // For second-order derivative:
  Eigen::MatrixXd all_Hessian_nu_n_c_equal_0 = Eigen::MatrixXd::Zero(max_h,1) ;
  Eigen::MatrixXd all_Hessian_nu_n_c         = Eigen::MatrixXd::Zero(max_h,1) ;
  
  Eigen::MatrixXd F_c_equal_0   = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd F_c           = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd Phi_aux_s_bar = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd phi_aux_s_bar = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd Phi_aux_s     = Eigen::MatrixXd::Zero(T,max_h) ;
  
  Eigen::MatrixXd cum_F_c_equal_0       = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd cum_F_c               = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd cum_F_c_equal_0_overN = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd cum_F_c_overN         = Eigen::MatrixXd::Zero(T,max_h) ;
  
  Eigen::MatrixXd Matrix_i_bar = Eigen::MatrixXd::Constant(T,max_h,i_bar) ;
  
  // Note: the following objects are relevant only if T == 1:
  // (used for computation of gradients and second-order derivatives)
  Eigen::MatrixXd cum_dF_c_equal_0       = Eigen::MatrixXd::Zero(max_h,n_x) ;
  Eigen::MatrixXd cum_dF_c               = Eigen::MatrixXd::Zero(max_h,n_x) ;
  Eigen::MatrixXd cum_dF_c_equal_0_overN = Eigen::MatrixXd::Zero(max_h,n_x) ;
  Eigen::MatrixXd cum_dF_c_overN         = Eigen::MatrixXd::Zero(max_h,n_x) ;
  Eigen::MatrixXd cum_d2F_c_equal_0       = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
  Eigen::MatrixXd cum_d2F_c               = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
  Eigen::MatrixXd cum_d2F_c_equal_0_overN = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
  Eigen::MatrixXd cum_d2F_c_overN         = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
  
  aux_an_bnX    = vec1T * all_a_n.transpose()     + X * all_b_n.transpose() - Matrix_i_bar ;
  aux_abarn_bnX = vec1T * all_a_bar_n.transpose() + X * all_b_n.transpose() - Matrix_i_bar ;
  aux_s         = aux_an_bnX.array()    * ((vec1T * all_sigma_n.transpose()).cwiseInverse()).array() ;
  aux_s_bar     = aux_abarn_bnX.array() * ((vec1T * all_sigma_n.transpose()).cwiseInverse()).array() ;
  
  aux_sigma = vec1T * all_sigma_n.transpose() ;
  aux_pi    = vec1T * all_c_dot_n.transpose() + X * all_c_n.transpose() ;
  
  F_c_equal_0 = Matrix_i_bar.array() + aux_sigma.array() * g(aux_s).array() ;
  
  Phi_aux_s_bar = pnorm_cpp(aux_s_bar) ;
  Phi_aux_s     = pnorm_cpp(aux_s) ;
  
  Eigen::MatrixXd auxaux    = Eigen::MatrixXd::Zero(max_h,1) ;
  auxaux = aux_SS.array()*(all_sigma_n.cwiseInverse()).array() ;
  phi_aux_s_bar = dnorm_cpp(aux_s_bar).array() * (vec1T * auxaux.transpose()).array() ;
  
  // For second-order derivative (IRF only):
  Eigen::MatrixXd phi_s_Hessian = Eigen::MatrixXd::Zero(max_h,1) ;
  phi_s_Hessian = dnorm_cpp(aux_s).array() * (all_sigma_n.cwiseInverse()).array() ;
  Eigen::MatrixXd phi_aux_s_bar_Hessian = Eigen::MatrixXd::Zero(max_h,1) ;
  phi_aux_s_bar_Hessian = dnorm_cpp(aux_s_bar).array() * (vec1T * auxaux.transpose()).array() * 
    (all_sigma_n.cwiseInverse()).array() * aux_s_bar.array() ;
  
  Eigen::MatrixXd vec_i_bar = Eigen::MatrixXd::Constant(T,1,i_bar) ;
  
  Eigen::MatrixXd indic_sigma0 = Eigen::MatrixXd::Zero(max_h,1) ;
  for (int i = 0; i < max_h; i++){
    if(all_sigma_n(i,0)==0){
      indic_sigma0(i) = 1 ;
      
      F_c_equal_0.col(i)   = vec_i_bar + pmax_cpp(aux_an_bnX.col(i), 0) ;
      // Phi_aux_s_bar.col(i) = positiveIndicator(aux_s_bar.col(i)) ;
      // Phi_aux_s.col(i)     = positiveIndicator(aux_s.col(i)) ;
      Phi_aux_s_bar.col(i) = vec1N ;
      Phi_aux_s.col(i)     = vec1N ;
      phi_aux_s_bar.col(i) = vec0N ;
      
      phi_s_Hessian(i,0)         = 0 ;
      phi_aux_s_bar_Hessian(i,0) = 0 ;
    }
  }
  
  // For real rates, add two additional terms:
  F_c = Phi_aux_s_bar.array() * (vec1T * aux_SS.transpose()).array() ;
  F_c = F_c + F_c_equal_0 - aux_pi ;
  
  // This is used to compute dF:
  all_nu_n_c_equal_0 = Phi_aux_s ;
  all_nu_n_c         = Phi_aux_s + phi_aux_s_bar ;
  
  // This is used to compute d2F (for IRFs):
  all_Hessian_nu_n_c_equal_0 = phi_s_Hessian ;
  all_Hessian_nu_n_c         = phi_s_Hessian - phi_aux_s_bar_Hessian ;
  
  // If X is a row vector, compute gradient of F_c's w.r.t. X:
  Eigen::MatrixXd dF_c         = Eigen::MatrixXd::Zero(max_h,n_x) ;
  Eigen::MatrixXd dF_c_equal_0 = Eigen::MatrixXd::Zero(max_h,n_x) ;
  // If X is a row vector, compute second-order derivative of F_c's w.r.t. X:
  // ()This is used only for IRFs)
  Eigen::MatrixXd d2F_c         = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
  Eigen::MatrixXd d2F_c_equal_0 = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
  
  if(T == 1){
    // in that case, the "all_nu_n" matrices are of dimension 1 x max_h (i.e., T == 1)
    Eigen::MatrixXd Diag_all_nu_n_c         = all_nu_n_c.asDiagonal() ;
    Eigen::MatrixXd Diag_all_nu_n_c_equal_0 = all_nu_n_c_equal_0.asDiagonal() ;
    
    dF_c_equal_0 = Diag_all_nu_n_c_equal_0 * all_b_n ;
    dF_c         = Diag_all_nu_n_c         * all_b_n - all_c_n ;
    
    //For second-order derivatives:
    Eigen::MatrixXd all_bb_n = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
    all_bb_n = kronecker_cpp(all_b_n,vec1x.transpose()).array() * 
      kronecker_cpp(vec1x.transpose(),all_b_n).array() ;
    
    Eigen::MatrixXd Diag_Hessian_all_nu_n_c         = all_Hessian_nu_n_c.asDiagonal() ;
    Eigen::MatrixXd Diag_Hessian_all_nu_n_c_equal_0 = all_Hessian_nu_n_c_equal_0.asDiagonal() ;
    
    d2F_c_equal_0 = Diag_Hessian_all_nu_n_c_equal_0 * all_bb_n ;
    d2F_c         = Diag_Hessian_all_nu_n_c         * all_bb_n ;
  }
  
  // Compute yields and gradients : --------------------------------------------
  cum_F_c_equal_0_overN = F_c_equal_0 ; // real yields
  cum_F_c_overN         = F_c ; // nominal yields
  cum_F_c_equal_0       = F_c_equal_0 ; // real yields
  cum_F_c               = F_c ; // nominal yields
  
  cum_dF_c_equal_0       = dF_c_equal_0 ;
  cum_dF_c               = dF_c ;
  cum_dF_c_equal_0_overN = dF_c_equal_0 ;
  cum_dF_c_overN         = dF_c ;
  
  cum_d2F_c_equal_0       = d2F_c_equal_0 ;
  cum_d2F_c               = d2F_c ;
  cum_d2F_c_equal_0_overN = d2F_c_equal_0 ;
  cum_d2F_c_overN         = d2F_c ;
  
  for (int i = 2; i <= max_h; i++){
    // Yields:
    cum_F_c_equal_0.col(i-1) = cum_F_c_equal_0.col(i-2) + F_c_equal_0.col(i-1) ;
    cum_F_c.col(i-1)         = cum_F_c.col(i-2)         + F_c.col(i-1) ;
    
    cum_F_c_equal_0_overN.col(i-1) = cum_F_c_equal_0.col(i-1)/i ;
    cum_F_c_overN.col(i-1)         = cum_F_c.col(i-1)/i ;
    
    // Yields gradient (meaningful only if T == 1):
    cum_dF_c_equal_0.row(i-1) = cum_dF_c_equal_0.row(i-2) + dF_c_equal_0.row(i-1) ;
    cum_dF_c.row(i-1)         = cum_dF_c.row(i-2)         + dF_c.row(i-1) ;
    
    cum_dF_c_equal_0_overN.row(i-1) = cum_dF_c_equal_0.row(i-1)/i ;
    cum_dF_c_overN.row(i-1)         = cum_dF_c.row(i-1)/i ;
    
    // Yields second-order derivatives (meaningful only if T == 1):
    cum_d2F_c_equal_0.row(i-1) = cum_d2F_c_equal_0.row(i-2) + d2F_c_equal_0.row(i-1) ;
    cum_d2F_c.row(i-1)         = cum_d2F_c.row(i-2)         + d2F_c.row(i-1) ;
    
    cum_d2F_c_equal_0_overN.row(i-1) = cum_d2F_c_equal_0.row(i-1)/i ;
    cum_d2F_c_overN.row(i-1)         = cum_d2F_c.row(i-1)/i ;
  }
  
  return List::create(Named("F_c_equal_0")            = F_c_equal_0,
                      Named("F_c")                    = F_c,
                      Named("all_b_n")                = all_b_n,
                      Named("all_c_n")                = all_c_n,
                      Named("all_nu_n_c")             = all_nu_n_c,
                      Named("all_nu_n_c_equal_0")     = all_nu_n_c_equal_0,
                      Named("cum_F_c_equal_0_overN")  = cum_F_c_equal_0_overN,
                      Named("cum_F_c_overN")          = cum_F_c_overN,
                      Named("cum_dF_c_equal_0_overN") = cum_dF_c_equal_0_overN,
                      Named("cum_dF_c")         = cum_dF_c,
                      Named("cum_dF_c_overN")         = cum_dF_c_overN,
                      Named("cum_d2F_c_equal_0_overN") = cum_d2F_c_equal_0_overN,
                      Named("cum_d2F_c_overN")         = cum_d2F_c_overN,
                      Named("aux_s") = aux_s) ;
}





// [[Rcpp::export]]
Rcpp::List apply_Lemma_cpp_from_RES(const Rcpp::List      & RES_aux,
                                    const Eigen::MatrixXd & X,
                                    const Eigen::MatrixXd & i_bar,
                                    const int max_h){
  
  Eigen::MatrixXd all_a_n       = RES_aux("all_a_n") ;
  Eigen::MatrixXd all_a_bar_n   = RES_aux("all_a_bar_n") ;
  Eigen::MatrixXd all_c_dot_n   = RES_aux("all_c_dot_n") ;
  Eigen::MatrixXd all_sigma_n   = RES_aux("all_sigma_n") ;
  Eigen::MatrixXd all_b_n       = RES_aux("all_b_n") ;
  Eigen::MatrixXd all_c_n       = RES_aux("all_c_n") ;
  Eigen::MatrixXd all_b_tilde_n = RES_aux("all_b_tilde_n") ;
  Eigen::MatrixXd all_c_tilde_n = RES_aux("all_c_tilde_n") ;
  Eigen::MatrixXd aux_SS        = RES_aux("aux_SS") ;
  Eigen::MatrixXd vecSS         = RES_aux("vecSS") ;
  
  int n_x = all_b_n.cols() ;
  int T   = X.rows() ;
  
  Eigen::MatrixXd vec1N   = Eigen::MatrixXd::Ones(max_h,1) ;
  Eigen::MatrixXd vec0N   = Eigen::MatrixXd::Zero(max_h,1) ;
  Eigen::MatrixXd vec1T   = Eigen::MatrixXd::Ones(T,1) ;
  Eigen::MatrixXd vec1x   = Eigen::MatrixXd::Ones(n_x,1) ;
  Eigen::MatrixXd vec1xx  = Eigen::MatrixXd::Ones(n_x*n_x,1) ;
  
  Eigen::MatrixXd aux_an_bnX    = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd aux_abarn_bnX = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd aux_s         = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd aux_s_bar     = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd aux_sigma     = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd aux_pi        = Eigen::MatrixXd::Zero(T,max_h) ;
  
  // For derivative:
  Eigen::MatrixXd all_nu_n_c_equal_0 = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd all_nu_n_c         = Eigen::MatrixXd::Zero(T,max_h) ;
  // For second-order derivative:
  Eigen::MatrixXd all_Hessian_nu_n_c_equal_0 = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd all_Hessian_nu_n_c         = Eigen::MatrixXd::Zero(T,max_h) ;
  
  Eigen::MatrixXd F_c_equal_0   = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd F_c           = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd Phi_aux_s_bar = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd phi_aux_s_bar = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd Phi_aux_s     = Eigen::MatrixXd::Zero(T,max_h) ;
  
  Eigen::MatrixXd cum_F_c_equal_0       = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd cum_F_c               = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd cum_F_c_equal_0_overN = Eigen::MatrixXd::Zero(T,max_h) ;
  Eigen::MatrixXd cum_F_c_overN         = Eigen::MatrixXd::Zero(T,max_h) ;
  
  Eigen::MatrixXd Matrix_i_bar = Eigen::MatrixXd::Zero(T,max_h) ;
  Matrix_i_bar = i_bar * vec1N.transpose() ;
  
  // Note: the following objects are relevant only if T == 1:
  // (used for computation of gradients and second-order derivatives)
  Eigen::MatrixXd cum_dF_c_equal_0       = Eigen::MatrixXd::Zero(max_h,n_x) ;
  Eigen::MatrixXd cum_dF_c               = Eigen::MatrixXd::Zero(max_h,n_x) ;
  Eigen::MatrixXd cum_dF_c_equal_0_overN = Eigen::MatrixXd::Zero(max_h,n_x) ;
  Eigen::MatrixXd cum_dF_c_overN         = Eigen::MatrixXd::Zero(max_h,n_x) ;
  Eigen::MatrixXd cum_d2F_c_equal_0       = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
  Eigen::MatrixXd cum_d2F_c               = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
  Eigen::MatrixXd cum_d2F_c_equal_0_overN = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
  Eigen::MatrixXd cum_d2F_c_overN         = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
  
  aux_an_bnX    = vec1T * all_a_n.transpose()     + X * all_b_n.transpose() - Matrix_i_bar ;
  aux_abarn_bnX = vec1T * all_a_bar_n.transpose() + X * all_b_n.transpose() - Matrix_i_bar ;
  aux_s         = aux_an_bnX.array()    * ((vec1T * all_sigma_n.transpose()).cwiseInverse()).array() ;
  aux_s_bar     = aux_abarn_bnX.array() * ((vec1T * all_sigma_n.transpose()).cwiseInverse()).array() ;
  
  aux_sigma = vec1T * all_sigma_n.transpose() ;
  aux_pi    = vec1T * all_c_dot_n.transpose() + X * all_c_n.transpose() ;
  
  F_c_equal_0 = Matrix_i_bar.array() + aux_sigma.array() * g(aux_s).array() ;
  
  Phi_aux_s_bar = pnorm_cpp(aux_s_bar) ;
  Phi_aux_s     = pnorm_cpp(aux_s) ;
  
  Eigen::MatrixXd auxaux    = Eigen::MatrixXd::Zero(max_h,1) ;
  auxaux = aux_SS.array()*(all_sigma_n.cwiseInverse()).array() ;
  phi_aux_s_bar = dnorm_cpp(aux_s_bar).array() * (vec1T * auxaux.transpose()).array() ;
  
  // For second-order derivative (IRF only):
  Eigen::MatrixXd phi_s_Hessian = Eigen::MatrixXd::Zero(T,max_h) ;
  phi_s_Hessian = dnorm_cpp(aux_s).array() * ((vec1T * all_sigma_n.transpose()).cwiseInverse()).array() ;
  Eigen::MatrixXd phi_aux_s_bar_Hessian = Eigen::MatrixXd::Zero(T,max_h) ;
  phi_aux_s_bar_Hessian = dnorm_cpp(aux_s_bar).array() * (vec1T * auxaux.transpose()).array() *
    ((vec1T * all_sigma_n.transpose()).cwiseInverse()).array() * aux_s_bar.array() ;
  
  Eigen::MatrixXd vec_i_bar = i_bar ;
  
  Eigen::MatrixXd indic_sigma0 = Eigen::MatrixXd::Zero(max_h,1) ;
  for (int i = 0; i < max_h; i++){
    if(all_sigma_n(i,0)==0){
      indic_sigma0(i) = 1 ;
      
      F_c_equal_0.col(i)   = vec_i_bar + pmax_cpp(aux_an_bnX.col(i), 0) ;
      // Phi_aux_s_bar.col(i) = positiveIndicator(aux_s_bar.col(i)) ;
      // Phi_aux_s.col(i)     = positiveIndicator(aux_s.col(i)) ;
      Phi_aux_s_bar.col(i) = vec1N ;
      Phi_aux_s.col(i)     = vec1N ;
      phi_aux_s_bar.col(i) = vec0N ;
      
      phi_s_Hessian(i,0)         = 0 ;
      phi_aux_s_bar_Hessian(i,0) = 0 ;
    }
  }
  
  // For real rates, add two additional terms:
  F_c = Phi_aux_s_bar.array() * (vec1T * aux_SS.transpose()).array() ;
  F_c = F_c + F_c_equal_0 - aux_pi ;
  
  // This is used to compute dF:
  all_nu_n_c_equal_0 = Phi_aux_s ;
  all_nu_n_c         = Phi_aux_s + phi_aux_s_bar ;
  
  // // This is used to compute d2F (for IRFs):
  all_Hessian_nu_n_c_equal_0 = phi_s_Hessian ;
  all_Hessian_nu_n_c         = phi_s_Hessian - phi_aux_s_bar_Hessian ;
  
  // If X is a row vector, compute gradient of F_c's w.r.t. X:
  Eigen::MatrixXd dF_c         = Eigen::MatrixXd::Zero(max_h,n_x) ;
  Eigen::MatrixXd dF_c_equal_0 = Eigen::MatrixXd::Zero(max_h,n_x) ;
  // If X is a row vector, compute second-order derivative of F_c's w.r.t. X:
  // ()This is used only for IRFs)
  Eigen::MatrixXd d2F_c         = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
  Eigen::MatrixXd d2F_c_equal_0 = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
  
  if(T == 1){
    // in that case, the "all_nu_n" matrices are of dimension 1 x max_h (i.e., T == 1)
    Eigen::MatrixXd Diag_all_nu_n_c         = all_nu_n_c.asDiagonal() ;
    Eigen::MatrixXd Diag_all_nu_n_c_equal_0 = all_nu_n_c_equal_0.asDiagonal() ;
    
    dF_c_equal_0 = Diag_all_nu_n_c_equal_0 * all_b_n ;
    dF_c         = Diag_all_nu_n_c         * all_b_n - all_c_n ;
    
    //For second-order derivatives:
    Eigen::MatrixXd all_bb_n = Eigen::MatrixXd::Zero(max_h,n_x*n_x) ;
    all_bb_n = kronecker_cpp(all_b_n,vec1x.transpose()).array() *
      kronecker_cpp(vec1x.transpose(),all_b_n).array() ;
    
    Eigen::MatrixXd Diag_Hessian_all_nu_n_c         = all_Hessian_nu_n_c.asDiagonal() ;
    Eigen::MatrixXd Diag_Hessian_all_nu_n_c_equal_0 = all_Hessian_nu_n_c_equal_0.asDiagonal() ;
    
    d2F_c_equal_0 = Diag_Hessian_all_nu_n_c_equal_0 * all_bb_n ;
    d2F_c         = Diag_Hessian_all_nu_n_c         * all_bb_n ;
  }
  
  // Compute yields and gradients : --------------------------------------------
  cum_F_c_equal_0_overN = F_c_equal_0 ; // real yields
  cum_F_c_overN         = F_c ; // nominal yields
  cum_F_c_equal_0       = F_c_equal_0 ; // real yields
  cum_F_c               = F_c ; // nominal yields
  
  cum_dF_c_equal_0       = dF_c_equal_0 ;
  cum_dF_c               = dF_c ;
  cum_dF_c_equal_0_overN = dF_c_equal_0 ;
  cum_dF_c_overN         = dF_c ;
  
  cum_d2F_c_equal_0       = d2F_c_equal_0 ;
  cum_d2F_c               = d2F_c ;
  cum_d2F_c_equal_0_overN = d2F_c_equal_0 ;
  cum_d2F_c_overN         = d2F_c ;
  
  for (int i = 2; i <= max_h; i++){
    // Yields:
    cum_F_c_equal_0.col(i-1) = cum_F_c_equal_0.col(i-2) + F_c_equal_0.col(i-1) ;
    cum_F_c.col(i-1)         = cum_F_c.col(i-2)         + F_c.col(i-1) ;
    
    cum_F_c_equal_0_overN.col(i-1) = cum_F_c_equal_0.col(i-1)/i ;
    cum_F_c_overN.col(i-1)         = cum_F_c.col(i-1)/i ;
    
    // Yields gradient (meaningful only if T == 1):
    cum_dF_c_equal_0.row(i-1) = cum_dF_c_equal_0.row(i-2) + dF_c_equal_0.row(i-1) ;
    cum_dF_c.row(i-1)         = cum_dF_c.row(i-2)         + dF_c.row(i-1) ;
    
    cum_dF_c_equal_0_overN.row(i-1) = cum_dF_c_equal_0.row(i-1)/i ;
    cum_dF_c_overN.row(i-1)         = cum_dF_c.row(i-1)/i ;
    
    // Yields second-order derivatives (meaningful only if T == 1):
    cum_d2F_c_equal_0.row(i-1) = cum_d2F_c_equal_0.row(i-2) + d2F_c_equal_0.row(i-1) ;
    cum_d2F_c.row(i-1)         = cum_d2F_c.row(i-2)         + d2F_c.row(i-1) ;
    
    cum_d2F_c_equal_0_overN.row(i-1) = cum_d2F_c_equal_0.row(i-1)/i ;
    cum_d2F_c_overN.row(i-1)         = cum_d2F_c.row(i-1)/i ;
  }
  
  return List::create(Named("F_c_equal_0")            = F_c_equal_0,
                      Named("F_c")                    = F_c,
                      Named("all_a_n")                = all_a_n,
                      Named("all_b_n")                = all_b_n,
                      Named("all_c_n")                = all_c_n,
                      Named("all_nu_n_c")             = all_nu_n_c,
                      Named("all_nu_n_c_equal_0")     = all_nu_n_c_equal_0,
                      Named("cum_F_c_equal_0_overN")  = cum_F_c_equal_0_overN,
                      Named("cum_F_c_overN")          = cum_F_c_overN,
                      Named("cum_dF_c_equal_0_overN") = cum_dF_c_equal_0_overN,
                      Named("cum_dF_c")         = cum_dF_c,
                      Named("cum_dF_c_overN")         = cum_dF_c_overN,
                      Named("cum_d2F_c_equal_0_overN") = cum_d2F_c_equal_0_overN,
                      Named("cum_d2F_c_overN")         = cum_d2F_c_overN) ;
}



// [[Rcpp::export]]
Rcpp::List apply_Lemma_cpp(const Eigen::MatrixXd & mu,
                           const Eigen::MatrixXd & Phi,
                           const Eigen::MatrixXd & Sigma,
                           const Eigen::MatrixXd & i_bar,
                           const double          & a,
                           const Eigen::MatrixXd & b,
                           const Eigen::MatrixXd & c,
                           const Eigen::MatrixXd & X,
                           const int max_h){
  
  Rcpp::List RES_aux = apply_Lemma_cpp_aux(mu,Phi,Sigma,a,b,c,max_h) ;
  
  Rcpp::List RES = apply_Lemma_cpp_from_RES(RES_aux,X,i_bar,max_h) ;
  
  return RES ;
}


// [[Rcpp::export]]
Rcpp::List prepare_forecasts(const Rcpp::List      & model_sol,
                             const Eigen::MatrixXd & horizon_fcsts
){
  
  Eigen::MatrixXd Phi    = model_sol("Phi") ;
  Eigen::MatrixXd mu     = model_sol("mu") ;
  Eigen::MatrixXd Sigma  = model_sol("Sigma") ;
  
  // Number of components in X:
  int n_x = Phi.rows() ;
  
  Eigen::MatrixXd PhiQ    = model_sol("PhiQ") ;
  Eigen::MatrixXd muQ     = model_sol("muQ") ;
  Eigen::MatrixXd SigmaSigma = Eigen::MatrixXd::Zero(n_x,n_x) ;
  SigmaSigma = Sigma * Sigma.transpose() ;
  
  
  // Inflation specification:
  double gamma0          = model_sol("gamma0") ;
  Eigen::MatrixXd gamma1 = model_sol("gamma1") ;
  
  // Real growth specification:
  double eta0          = model_sol("eta0") ;
  Eigen::MatrixXd eta1 = model_sol("eta1") ;
  
  // Shadow rate specification:
  double delta0          = model_sol("delta0") ;
  Eigen::MatrixXd delta1 = model_sol("delta1") ;
  
  // Compute specification of 10-year yield (linear specif, no max() constraint)
  int maturity_10y = 120 ;
  double A10 = - delta0 ;
  Eigen::MatrixXd B10 = Eigen::MatrixXd::Zero(n_x,1) ;
  B10 = - delta1 ;
  Eigen::MatrixXd aux = Eigen::MatrixXd::Zero(n_x,1) ;
  double a10 ;
  Eigen::MatrixXd b10 = Eigen::MatrixXd::Zero(n_x,1) ;
  for(int h = 1; h < maturity_10y; h++){
    aux = B10 - delta1 ;
    A10 = - delta0 + A10 + (aux.transpose()*muQ)(0,0)  + 
      0.5 * (aux.transpose()*SigmaSigma*aux)(0,0) ;
    B10 = (PhiQ.transpose()) * aux ;
  }
  a10 = - A10 / maturity_10y ;
  Eigen::MatrixXd vec_10y = Eigen::MatrixXd::Constant(1,n_x,maturity_10y) ;
  Eigen::MatrixXd vec_10y_1 = vec_10y.cwiseInverse() ;
  b10.col(0) = - B10.array() * vec_10y_1.array() ;
  
  // Rcout << "B10 : " << B10 << "\n";
  // Rcout << "A10 : " << A10 << "\n";
  // Rcout << "b10 : " << b10 << "\n";
  // Rcout << "a10 : " << a10 << "\n";
  
  // Computation forecasts for all horizons: -----------------------------------
  int nb_horiz_fcsts    = horizon_fcsts.rows() ;
  int max_horizon_fcsts = horizon_fcsts(nb_horiz_fcsts - 1) ;
  
  Eigen::MatrixXd fcst_SR_csts      = Eigen::MatrixXd::Zero(max_horizon_fcsts,1) ;
  Eigen::MatrixXd fcst_SR_mult      = Eigen::MatrixXd::Zero(max_horizon_fcsts,n_x) ;
  Eigen::MatrixXd fcst_SR_sigma     = Eigen::MatrixXd::Zero(max_horizon_fcsts,1) ;
  
  Eigen::MatrixXd fcst_avg_Infl_csts    = Eigen::MatrixXd::Zero(max_horizon_fcsts,1) ;
  Eigen::MatrixXd fcst_avg_Infl_mult    = Eigen::MatrixXd::Zero(max_horizon_fcsts,n_x) ;
  Eigen::MatrixXd fcst_avg_Growth_csts  = Eigen::MatrixXd::Zero(max_horizon_fcsts,1) ;
  Eigen::MatrixXd fcst_avg_Growth_mult  = Eigen::MatrixXd::Zero(max_horizon_fcsts,n_x) ;
  Eigen::MatrixXd fcst_avg_Y10_csts     = Eigen::MatrixXd::Zero(max_horizon_fcsts,1) ;
  Eigen::MatrixXd fcst_avg_Y10_mult     = Eigen::MatrixXd::Zero(max_horizon_fcsts,n_x) ;
  
  Eigen::MatrixXd mu_h     = Eigen::MatrixXd::Zero(n_x,1) ;
  Eigen::MatrixXd Phi_h    = Eigen::MatrixXd::Identity(n_x,n_x) ;
  Eigen::MatrixXd Var_h    = Eigen::MatrixXd::Zero(n_x,n_x) ;
  Eigen::MatrixXd EXh      = Eigen::MatrixXd::Zero(n_x,1) ;
  Eigen::MatrixXd EcumXh   = Eigen::MatrixXd::Zero(n_x,1) ;
  Eigen::MatrixXd cummu_h  = Eigen::MatrixXd::Zero(n_x,1) ;
  Eigen::MatrixXd cumPhi_h = Eigen::MatrixXd::Zero(n_x,n_x) ;
  
  Eigen::MatrixXd vec_h   = Eigen::MatrixXd::Zero(1,n_x) ;
  Eigen::MatrixXd vec_h_1 = Eigen::MatrixXd::Zero(1,n_x) ;
  
  for(int h = 1; h <= max_horizon_fcsts; h++){
    Var_h    = Var_h + Phi_h * Sigma * Sigma.transpose() * Phi_h.transpose() ;
    mu_h     = mu_h + Phi_h * mu ;
    Phi_h    = Phi_h * Phi ;
    cummu_h  = cummu_h  + mu_h ;
    cumPhi_h = cumPhi_h + Phi_h ;
    
    fcst_SR_csts(h-1,0)   = delta0 + (delta1.transpose() * mu_h)(0,0) ;
    fcst_SR_mult.row(h-1) = delta1.transpose() * Phi_h ;
    fcst_SR_sigma(h-1,0)  = sqrt(((kronecker_cpp(delta1,delta1).transpose()) * Vec(Var_h))(0,0)) ;
    
    vec_h = Eigen::MatrixXd::Constant(1,n_x,h) ;
    vec_h_1 = vec_h.cwiseInverse() ;
    
    fcst_avg_Growth_csts(h-1,0) = eta0   + ((eta1.transpose()   * cummu_h)(0,0))/h ;
    fcst_avg_Infl_csts(h-1,0)   = gamma0 + ((gamma1.transpose() * cummu_h)(0,0))/h ;
    fcst_avg_Y10_csts(h-1,0)    = a10    + ((b10.transpose()    * cummu_h)(0,0))/h ;
    
    fcst_avg_Growth_mult.row(h-1) = (eta1.transpose()   * cumPhi_h).array() * vec_h_1.array();
    fcst_avg_Infl_mult.row(h-1)   = (gamma1.transpose() * cumPhi_h).array() * vec_h_1.array();
    fcst_avg_Y10_mult.row(h-1)    = (b10.transpose()    * cumPhi_h).array() * vec_h_1.array();
  }
  
  // Rcout << "fcst_avg_Y10_mult : " << fcst_avg_Y10_mult << "\n";
  // Rcout << "fcst_avg_Y10_csts : " << fcst_avg_Y10_csts << "\n";
  
  
  return List::create(Named("fcst_SR_csts")         = fcst_SR_csts,
                      Named("fcst_SR_mult")         = fcst_SR_mult,
                      Named("fcst_SR_sigma")        = fcst_SR_sigma,
                      Named("fcst_avg_Infl_csts")   = fcst_avg_Infl_csts,
                      Named("fcst_avg_Infl_mult")   = fcst_avg_Infl_mult,
                      Named("fcst_avg_Growth_csts") = fcst_avg_Growth_csts,
                      Named("fcst_avg_Growth_mult") = fcst_avg_Growth_mult,
                      Named("fcst_avg_Y10_csts")    = fcst_avg_Y10_csts,
                      Named("fcst_avg_Y10_mult")    = fcst_avg_Y10_mult) ;
}




// [[Rcpp::export]]
Rcpp::List f_measurement(const  Rcpp::List      & model_sol,
                         const  Rcpp::List      & RES_aux,
                         const  Rcpp::List      & RES_fcsts,
                         const  Eigen::MatrixXd & maturity_nomi_yds,
                         const  Eigen::MatrixXd & maturity_real_yds,
                         const  Eigen::MatrixXd & horizon_fcsts,
                         const  Eigen::MatrixXd & X,
                         const  Eigen::MatrixXd & i_bar
){
  // Note: X has to be a matrix with only one row.
  // i_bar has to be a matrix (1,1)!!!!!
  
  //Rcout << "Obs_pred : " << X << "\n";
  
  // RES_aux is the output of function "apply_Lemma_cpp_aux"
  int nb_mat_nomi_yds = maturity_nomi_yds.rows() ;
  int nb_mat_real_yds = maturity_real_yds.rows() ;
  double max_h        = maturity_nomi_yds(nb_mat_nomi_yds-1,0) ;
  
  //double i_bar = model_sol("i_bar") ;
  
  Rcpp::List RES = apply_Lemma_cpp_from_RES(RES_aux,X,i_bar,max_h) ;
  //Rcpp::List RES = RES_aux ;
  
  Eigen::MatrixXd Phi    = model_sol("Phi") ;
  Eigen::MatrixXd mu     = model_sol("mu") ;
  Eigen::MatrixXd Sigma  = model_sol("Sigma") ;
  
  Eigen::MatrixXd fcst_SR_csts         = RES_fcsts("fcst_SR_csts") ;
  Eigen::MatrixXd fcst_SR_mult         = RES_fcsts("fcst_SR_mult") ;
  Eigen::MatrixXd fcst_SR_sigma        = RES_fcsts("fcst_SR_sigma") ;
  Eigen::MatrixXd fcst_avg_Infl_csts   = RES_fcsts("fcst_avg_Infl_csts") ;
  Eigen::MatrixXd fcst_avg_Infl_mult   = RES_fcsts("fcst_avg_Infl_mult") ;
  Eigen::MatrixXd fcst_avg_Growth_csts = RES_fcsts("fcst_avg_Growth_csts") ;
  Eigen::MatrixXd fcst_avg_Growth_mult = RES_fcsts("fcst_avg_Growth_mult") ;
  Eigen::MatrixXd fcst_avg_Y10_csts    = RES_fcsts("fcst_avg_Y10_csts") ;
  Eigen::MatrixXd fcst_avg_Y10_mult    = RES_fcsts("fcst_avg_Y10_mult") ;
  
  // Number of components in X:
  int n_x = Phi.rows() ;
  
  // double delta  = model_sol("delta") ;
  
  // Inflation specification:
  double gamma0          = model_sol("gamma0") ;
  Eigen::MatrixXd gamma1 = model_sol("gamma1") ;
  
  // Real growth specification:
  double eta0          = model_sol("eta0") ;
  Eigen::MatrixXd eta1 = model_sol("eta1") ;
  
  // Shadow rate specification:
  double delta0          = model_sol("delta0") ;
  Eigen::MatrixXd delta1 = model_sol("delta1") ;
  
  // Output gap specification:
  Eigen::MatrixXd eta1z = model_sol("eta1z") ;
  
  // PTR specification:
  Eigen::MatrixXd eta1pistar = model_sol("eta1pistar") ;
  
  Eigen::MatrixXd seq_max_h = Eigen::MatrixXd::Zero(max_h,1) ;
  seq_max_h = seqEigen(1,max_h) ;
  
  // Extract bond yields: ------------------------------------------------------
  Eigen::MatrixXd yds_nomi = Eigen::MatrixXd::Zero(1,max_h) ;
  Eigen::MatrixXd yds_real = Eigen::MatrixXd::Zero(1,max_h) ;
  yds_real = RES("cum_F_c_overN") ;
  yds_nomi = RES("cum_F_c_equal_0_overN") ;
  
  Eigen::MatrixXd dyds_nomi = Eigen::MatrixXd::Zero(max_h,n_x) ;
  Eigen::MatrixXd dyds_real = Eigen::MatrixXd::Zero(max_h,n_x) ;
  dyds_real = RES("cum_dF_c_overN") ;
  dyds_nomi = RES("cum_dF_c_equal_0_overN") ;
  
  // Computation forecasts for all horizons: -----------------------------------
  int nb_horiz_fcsts  = horizon_fcsts.rows() ;
  int max_horiz_fcsts = horizon_fcsts(nb_horiz_fcsts-1,0) ;
  
  double E_st     = 0;
  double sigma_st = 0;
  double E_it     = 0;
  double E_cumit  = 0;
  
  Eigen::MatrixXd E_avgSTRt  = Eigen::MatrixXd::Zero(max_horiz_fcsts,1) ;
  Eigen::MatrixXd dE_avgSTRt = Eigen::MatrixXd::Zero(max_horiz_fcsts,n_x) ;
  Eigen::MatrixXd dE_STRt    = Eigen::MatrixXd::Zero(1,n_x) ;
  Eigen::MatrixXd dE_cumSTRt = Eigen::MatrixXd::Zero(1,n_x) ;
  
  for(int h = 1; h <= max_horiz_fcsts; h++){
    // Short Term Rate:
    E_st   = - i_bar(0,0) + fcst_SR_csts(h-1,0) + (fcst_SR_mult.row(h-1) * X.transpose())(0,0) ;
    sigma_st = fcst_SR_sigma(h-1,0) ;
    double auxaux = E_st/sigma_st ;
    double Phi_EoverS = pnorm_scalar_cpp(auxaux) ;
    double phi_EoverS = dnorm_scalar_cpp(auxaux) ;
    Eigen::MatrixXd vecPhi_EoverS = Eigen::MatrixXd::Constant(1,n_x,Phi_EoverS) ;
    E_it    = i_bar(0,0) + Phi_EoverS * E_st + phi_EoverS * sigma_st ;
    E_cumit = E_cumit + E_it ;
    E_avgSTRt(h-1,0) = E_cumit / h ;
    
    dE_STRt    = vecPhi_EoverS.array() * (fcst_SR_mult.row(h-1)).array() ;
    dE_cumSTRt = dE_cumSTRt + dE_STRt ;
    Eigen::MatrixXd vec_h = Eigen::MatrixXd::Constant(1,n_x,h) ;
    dE_avgSTRt.row(h-1) = dE_cumSTRt.array() * (vec_h.cwiseInverse()).array() ;
  }
  
  // Point values:
  Eigen::MatrixXd infl_fcst   = Eigen::MatrixXd::Zero(nb_horiz_fcsts,1) ;
  Eigen::MatrixXd growth_fcst = Eigen::MatrixXd::Zero(nb_horiz_fcsts,1) ;
  Eigen::MatrixXd str_fcst    = Eigen::MatrixXd::Zero(nb_horiz_fcsts,1) ;
  Eigen::MatrixXd y10_fcst    = Eigen::MatrixXd::Zero(nb_horiz_fcsts,1) ;
  // Gradients:
  Eigen::MatrixXd dinfl_fcst   = Eigen::MatrixXd::Zero(nb_horiz_fcsts,n_x) ;
  Eigen::MatrixXd dgrowth_fcst = Eigen::MatrixXd::Zero(nb_horiz_fcsts,n_x) ;
  Eigen::MatrixXd dstr_fcst    = Eigen::MatrixXd::Zero(nb_horiz_fcsts,n_x) ;
  Eigen::MatrixXd dy10_fcst    = Eigen::MatrixXd::Zero(nb_horiz_fcsts,n_x) ;
  
  // Inflation, growth, STR, and Y10 forecasts for appropriate horizons:
  for(int nbh = 1; nbh <= nb_horiz_fcsts; nbh++){
    int h = horizon_fcsts(nbh - 1,0) ;
    infl_fcst(nbh-1,0)   = fcst_avg_Infl_csts(h-1,0)   + (fcst_avg_Infl_mult.row(h-1)   * X.transpose())(0,0) ;
    growth_fcst(nbh-1,0) = fcst_avg_Growth_csts(h-1,0) + (fcst_avg_Growth_mult.row(h-1) * X.transpose())(0,0) ;
    y10_fcst(nbh-1,0)    = fcst_avg_Y10_csts(h-1,0)    + (fcst_avg_Y10_mult.row(h-1) * X.transpose())(0,0) ;
    str_fcst(nbh-1,0)    = E_avgSTRt(h-1,0) ;

    dinfl_fcst.row(nbh-1)   = fcst_avg_Infl_mult.row(h-1) ;
    dgrowth_fcst.row(nbh-1) = fcst_avg_Growth_mult.row(h-1) ;
    dy10_fcst.row(nbh-1)    = fcst_avg_Y10_mult.row(h-1) ;
    dstr_fcst.row(nbh-1)    = dE_avgSTRt.row(h-1) ;
    
    // Rcout << "h : " << h << "\n";
    // Rcout << "fcst_avg_Infl_csts : " << fcst_avg_Infl_csts(h-1,0) << "\n";
    // Rcout << "fcst_avg_Growth_csts : " << fcst_avg_Growth_csts(h-1,0) << "\n";
    
  }
  
  int nb_obs_macro = 5 ; // real growth, inflation, short-term rate, output gap, PTR
  int nb_obs       = nb_obs_macro + nb_mat_nomi_yds + nb_mat_real_yds + 4*nb_horiz_fcsts ;
  // Note regarding above: 4*nb_horiz_fcsts because 4 variables: (infl & growth & TBILL & Y10)
  
  Eigen::MatrixXd fX  = Eigen::MatrixXd::Zero(nb_obs,1) ; 
  Eigen::MatrixXd dfX = Eigen::MatrixXd::Zero(nb_obs,n_x) ;
  
  fX(0,0)      = eta0 + (X * eta1)(0,0) ;
  dfX.row(0) = eta1.transpose() ;
  fX(1,0)      = gamma0 + (X * gamma1)(0,0) ;
  dfX.row(1) = gamma1.transpose() ;
  fX(2,0)      = max_cpp(delta0 + (X * delta1)(0,0),i_bar(0,0)) ;
  if(delta0 + (X * delta1)(0,0)>i_bar(0,0)){
    dfX.row(2) = delta1.transpose() ;
  }
  fX(3,0)      = (X * eta1z)(0,0) ;
  dfX.row(3) = eta1z.transpose() ;
  fX(4,0)      = (X * eta1pistar)(0,0) ;
  dfX.row(4) = eta1pistar.transpose() ;
  
  int counts ;
  // Nominal yields:
  counts = nb_obs_macro ;
  int indic_yd ;
  for(int i = 1; i <= nb_mat_nomi_yds; i++){
    indic_yd = maturity_nomi_yds(i-1,0) ;
    fX(counts + i - 1,0)    = yds_nomi(0,indic_yd-1) ;
    dfX.row(counts + i - 1) = dyds_nomi.row(indic_yd-1) ;
  }
  // Real yields:
  counts = nb_obs_macro + nb_mat_nomi_yds ;
  for(int i = 1; i <= nb_mat_real_yds; i++){
    indic_yd = maturity_real_yds(i-1,0) ;
    fX(counts + i - 1,0)    = yds_real(0,indic_yd-1) ;
    dfX.row(counts + i - 1) = dyds_real.row(indic_yd-1) ;
  }
  // Growth forecasts:
  counts = nb_obs_macro + nb_mat_nomi_yds + nb_mat_real_yds ;
  for(int i = 1; i <= nb_horiz_fcsts; i++){
    fX(counts + i - 1,0)    = growth_fcst(i-1,0) ;
    dfX.row(counts + i - 1) = dgrowth_fcst.row(i-1) ;
  }
  // Inflation forecasts:
  counts = nb_obs_macro + nb_mat_nomi_yds + nb_mat_real_yds + nb_horiz_fcsts;
  for(int i = 1; i <= nb_horiz_fcsts; i++){
    fX(counts + i - 1,0)    = infl_fcst(i-1,0) ;
    dfX.row(counts + i - 1) = dinfl_fcst.row(i-1) ;
  }
  // Interest rate forecasts:
  counts = nb_obs_macro + nb_mat_nomi_yds + nb_mat_real_yds + 2*nb_horiz_fcsts;
  for(int i = 1; i <= nb_horiz_fcsts; i++){
    fX(counts + i - 1,0)    = str_fcst(i-1,0) ;
    dfX.row(counts + i - 1) = dstr_fcst.row(i-1) ;
  }
  // 10-year yield forecasts:
  counts = nb_obs_macro + nb_mat_nomi_yds + nb_mat_real_yds + 3*nb_horiz_fcsts;
  for(int i = 1; i <= nb_horiz_fcsts; i++){
    fX(counts + i - 1,0)    = y10_fcst(i-1,0) ;
    dfX.row(counts + i - 1) = dy10_fcst.row(i-1) ;
  }
  
  return List::create(Named("fX")  = fX,
                      Named("dfX") = dfX) ;
}



// [[Rcpp::export]]
Rcpp::List EKF_filter_cpp(const Rcpp::List & model_sol,
                          const Rcpp::List & list_stdv,
                          const Rcpp::List & DATA,
                          const int indic_compute_fitted){
  
  const Eigen::MatrixXd stdv_yds_nomi     = list_stdv("stdv_yds_nomi") ;
  const Eigen::MatrixXd stdv_yds_real     = list_stdv("stdv_yds_real") ;
  const Eigen::MatrixXd stdv_str_fcsts    = list_stdv("stdv_str_fcsts") ;
  const Eigen::MatrixXd stdv_infl_fcsts   = list_stdv("stdv_infl_fcsts") ;
  const Eigen::MatrixXd stdv_growth_fcsts = list_stdv("stdv_growth_fcsts") ;
  const Eigen::MatrixXd stdv_y10_fcsts    = list_stdv("stdv_y10_fcsts") ;
  const Eigen::MatrixXd stdv_macro        = list_stdv("stdv_macro") ;
  
  Eigen::MatrixXd maturity_nomi_yds  = DATA("maturity_nomi_yds") ;
  Eigen::MatrixXd maturity_real_yds  = DATA("maturity_real_yds") ;
  Eigen::MatrixXd horizon_fcsts      = DATA("horizon_fcsts") ;
  
  int nb_yds_nomi   = maturity_nomi_yds.rows() ;
  int nb_yds_real   = maturity_real_yds.rows() ;
  int nb_fcsts = horizon_fcsts.rows() ; // per forecasted variable
  int max_h    = maturity_nomi_yds(nb_yds_nomi - 1) ;
  
  Eigen::MatrixXd yds_nomi     = DATA("yields_nomi") ;
  Eigen::MatrixXd yds_real     = DATA("yields_real") ;
  Eigen::MatrixXd infl_fcsts   = DATA("infl_fcsts") ;
  Eigen::MatrixXd growth_fcsts = DATA("growth_fcsts") ;
  Eigen::MatrixXd str_fcsts    = DATA("str_fcsts") ;
  Eigen::MatrixXd y10_fcsts    = DATA("y10_fcsts") ;

  int time = yds_nomi.rows() ;
  int nb_obs_macro = 5 ; // 5 macro variables: d.x, d.y, pi, s, PTR
  int nb_obs  = nb_obs_macro + nb_yds_nomi + nb_yds_real + 4*nb_fcsts ; //4 fcsted variables
  
  Eigen::MatrixXd dy     = DATA("dy") ;
  Eigen::MatrixXd infl   = DATA("inflation") ;
  Eigen::MatrixXd str    = DATA("str") ;
  Eigen::MatrixXd z      = DATA("z") ;
  Eigen::MatrixXd PTR    = DATA("PTR") ;
  
  Eigen::MatrixXd all_i_bar = Eigen::MatrixXd::Zero(time, 1) ;
  all_i_bar = DATA("i_bar") ;
  Eigen::MatrixXd i_bar = Eigen::MatrixXd::Zero(1, 1) ;

  // Data that will be fitted by filter:
  Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(time, nb_obs) ;
  
  Y.col(0) = dy ;
  Y.col(1) = infl ;
  Y.col(2) = str ;
  Y.col(3) = z ;
  Y.col(4) = PTR ;
  
  Y.block(0,nb_obs_macro,time,nb_yds_nomi)             = yds_nomi ;
  Y.block(0,nb_obs_macro+nb_yds_nomi,time,nb_yds_real) = yds_real ;
  Y.block(0,nb_obs_macro+nb_yds_nomi+nb_yds_real,time,nb_fcsts)            = growth_fcsts ;
  Y.block(0,nb_obs_macro+nb_yds_nomi+nb_yds_real+nb_fcsts,time,nb_fcsts)   = infl_fcsts ;
  Y.block(0,nb_obs_macro+nb_yds_nomi+nb_yds_real+2*nb_fcsts,time,nb_fcsts) = str_fcsts ;
  Y.block(0,nb_obs_macro+nb_yds_nomi+nb_yds_real+3*nb_fcsts,time,nb_fcsts) = y10_fcsts ;
  
  Eigen::MatrixXd Sigma = model_sol("Sigma") ;
  int n_x   = Sigma.rows() ;
  int n_eps = Sigma.cols() ;
  Eigen::MatrixXd mu    = Eigen::MatrixXd::Zero(n_x, 1) ;
  Eigen::MatrixXd Phi   = Eigen::MatrixXd::Zero(n_x, n_x) ;
  Eigen::MatrixXd muQ   = Eigen::MatrixXd::Zero(n_x, 1) ;
  Eigen::MatrixXd PhiQ  = Eigen::MatrixXd::Zero(n_x, n_x) ;
  
  mu    = model_sol("mu") ;
  Phi   = model_sol("Phi") ;
  Sigma = model_sol("Sigma") ;
  muQ   = model_sol("muQ") ;
  PhiQ  = model_sol("PhiQ") ;
  
  // Bond pricing:
  
  double a = model_sol("delta0") ;
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(n_x, 1) ;
  Eigen::MatrixXd c = Eigen::MatrixXd::Zero(n_x, 1) ;
  b = model_sol("delta1_tilde") ;
  c = model_sol("gamma1") ;
  
  Rcpp::List RES_aux = apply_Lemma_cpp_aux(muQ,PhiQ,Sigma,
                                           a,b,c,max_h) ;
  
  Rcpp::List RES_fcsts = prepare_forecasts(model_sol,horizon_fcsts) ;
  
  Eigen::MatrixXd H     = Eigen::MatrixXd::Zero(n_x, n_eps) ;
  Eigen::MatrixXd N     = Eigen::MatrixXd::Zero(n_x, n_eps) ;
  Eigen::MatrixXd G     = Eigen::MatrixXd::Zero(nb_obs, n_x) ;
  Eigen::MatrixXd M     = Eigen::MatrixXd::Zero(nb_obs, nb_obs) ;
  Eigen::MatrixXd diagM                  = Eigen::MatrixXd::Zero(nb_obs, 1) ;
  Eigen::MatrixXd diag_stdv_macro        = Eigen::MatrixXd::Zero(nb_obs_macro, 1) ;
  Eigen::MatrixXd diag_stdv_yds_nomi     = Eigen::MatrixXd::Zero(nb_yds_nomi, 1) ;
  Eigen::MatrixXd diag_stdv_yds_real     = Eigen::MatrixXd::Zero(nb_yds_real, 1) ;
  Eigen::MatrixXd diag_stdv_growth_fcsts = Eigen::MatrixXd::Zero(nb_fcsts, 1) ;
  Eigen::MatrixXd diag_stdv_infl_fcsts   = Eigen::MatrixXd::Zero(nb_fcsts, 1) ;
  Eigen::MatrixXd diag_stdv_str_fcsts    = Eigen::MatrixXd::Zero(nb_fcsts, 1) ;
  Eigen::MatrixXd diag_stdv_y10_fcsts    = Eigen::MatrixXd::Zero(nb_fcsts, 1) ;
  
  int counts ;
  
  H = model_sol("Phi") ;
  N = Sigma ;
  
  // // Prepare M matrix
  // Eigen::MatrixXd Diag_stdv_macro    = stdv_macro.asDiagonal() ;
  // Eigen::MatrixXd diag_stdv_yds_nomi = Eigen::MatrixXd::Constant(nb_yds_nomi,1, stdv_yds_nomi) ;
  // Eigen::MatrixXd Diag_stdv_yds_nomi = diag_stdv_yds_nomi.asDiagonal() ;
  // Eigen::MatrixXd diag_stdv_yds_real = Eigen::MatrixXd::Constant(nb_yds_real,1, stdv_yds_real) ;
  // Eigen::MatrixXd Diag_stdv_yds_real = diag_stdv_yds_real.asDiagonal() ;
  // Eigen::MatrixXd Diag_stdv_growth_fcsts = stdv_growth_fcsts.asDiagonal() ;
  // Eigen::MatrixXd Diag_stdv_infl_fcsts   = stdv_infl_fcsts.asDiagonal() ;
  // Eigen::MatrixXd Diag_stdv_str_fcsts    = stdv_str_fcsts.asDiagonal() ;
  // M.block(0,0,nb_obs_macro,nb_obs_macro) = Diag_stdv_macro ;
  // int counts = nb_obs_macro;
  // M.block(counts,counts,nb_yds_nomi,nb_yds_nomi) = Diag_stdv_yds_nomi ;
  // counts = nb_obs_macro + nb_yds_nomi ;
  // M.block(counts,counts,nb_yds_real,nb_yds_real) = Diag_stdv_yds_real ;
  // counts = nb_obs_macro + nb_yds_nomi + nb_yds_real ;
  // M.block(counts,counts,nb_fcsts,nb_fcsts) = Diag_stdv_growth_fcsts ;
  // counts = nb_obs_macro + nb_yds_nomi + nb_yds_real + nb_fcsts ;
  // M.block(counts,counts,nb_fcsts,nb_fcsts) = Diag_stdv_infl_fcsts ;
  // counts = nb_obs_macro + nb_yds_nomi + nb_yds_real + 2*nb_fcsts ;
  // M.block(counts,counts,nb_fcsts,nb_fcsts) = Diag_stdv_str_fcsts ;
  
  Eigen::MatrixXd P0 = model_sol("VX") ;
  Eigen::MatrixXd W0 = model_sol("EX") ;
  
  // Building objects to retrieve
  //-----------------------------
  double pi = 3.141592653589793;
  
  Eigen::MatrixXd W_t_t          = Eigen::MatrixXd::Zero(n_x, time);
  Eigen::MatrixXd W_t_t_nomodif  = Eigen::MatrixXd::Zero(n_x, time);
  Eigen::MatrixXd Obs_fitted     = Eigen::MatrixXd::Zero(nb_obs, time);
  
  Rcpp::List P_t_t   = Rcpp::List(time);
  
  Eigen::MatrixXd loglik = Eigen::MatrixXd::Zero(time,1);
  
  // Initialize the Filter
  //----------------------
  Eigen::MatrixXd W = W0;
  Eigen::MatrixXd P = P0;
  
  // Rcout << "W0 : " << W << "\n";
  // Rcout << "P0 : " << P << "\n";
  
  Eigen::MatrixXd W_pred          = Eigen::MatrixXd::Zero(n_x,1);
  Eigen::MatrixXd Obs_pred        = Eigen::MatrixXd::Zero(nb_obs,1);
  Eigen::MatrixXd Obs_fitted_t    = Eigen::MatrixXd::Zero(nb_obs,1);
  
  Eigen::MatrixXd W_up  = Eigen::MatrixXd::Zero(n_x,1);
  Eigen::MatrixXd P_up  = Eigen::MatrixXd::Zero(n_x, n_x);
  
  Eigen::MatrixXd P_pred          = Eigen::MatrixXd::Zero(n_x, n_x);
  Eigen::MatrixXd M_pred          = Eigen::MatrixXd::Zero(nb_obs, nb_obs);
  
  Eigen::MatrixXd Condi_P_pred    = Eigen::MatrixXd::Zero(n_x, n_x);
  
  Eigen::MatrixXd X_t = Eigen::MatrixXd::Zero(1,n_x);
  
  Eigen::MatrixXd nut       = Eigen::MatrixXd::Zero(n_x,1);
  nut = model_sol("mu") ;
  Eigen::MatrixXd Yt        = Eigen::MatrixXd::Zero(nb_obs,1);
  Eigen::MatrixXd residual  = Eigen::MatrixXd::Zero(nb_obs,1);
  
  Rcpp::List list_Xt = Rcpp::List(time);
  Rcpp::List list_G  = Rcpp::List(time);
  
  Eigen::MatrixXd Gain ;
  
  // LAUNCH THE RECURSIONS
  //======================
  for(int t = 0; t < time; t++){
    //for(int t = 0; t < 2; t++){
    
    // Lower bound:
    i_bar(0,0) = all_i_bar(t,0) ;
    
    //nut  = (nu_t.row(t)).transpose() ;
    Yt   = (Y.row(t)).transpose() ;
    
    // PREDICTION STEP
    //----------------
    // Mean of latent factors
    W_pred = nut + H * W;
    
    // Predict the variance-covariance matrix
    P_pred = H * P * H.transpose() + N * N.transpose();
    
    X_t = W_pred.transpose() ;
    
    // if(t == 213){
    //   Rcout << "X : " << X_t << "\n";
    // }
    
    Rcpp::List RES_measurement = f_measurement(model_sol,
                                               RES_aux,
                                               RES_fcsts,
                                               maturity_nomi_yds,
                                               maturity_real_yds,
                                               horizon_fcsts,
                                               X_t,
                                               i_bar) ;
    G         = RES_measurement("dfX") ;
    list_G(t) = G ;
    
    // Construct M:
    diag_stdv_macro        = stdv_macro.row(t) ;
    diag_stdv_yds_nomi     = stdv_yds_nomi.row(t) ;
    diag_stdv_yds_real     = stdv_yds_real.row(t) ;
    diag_stdv_growth_fcsts = stdv_growth_fcsts.row(t) ;
    diag_stdv_infl_fcsts   = stdv_infl_fcsts.row(t) ;
    diag_stdv_str_fcsts    = stdv_str_fcsts.row(t) ;
    diag_stdv_y10_fcsts    = stdv_y10_fcsts.row(t) ;
    counts = 0 ;
    diagM.block(counts,0,nb_obs_macro,1) = diag_stdv_macro ;
    counts = counts + nb_obs_macro ;
    diagM.block(counts,0,nb_yds_nomi,1)  = diag_stdv_yds_nomi ;
    counts = counts + nb_yds_nomi ;
    diagM.block(counts,0,nb_yds_real,1)  = diag_stdv_yds_real ;
    counts = counts + nb_yds_real ;
    diagM.block(counts,0,nb_fcsts,1)     = diag_stdv_growth_fcsts ;
    counts = counts + nb_fcsts ;
    diagM.block(counts,0,nb_fcsts,1)     = diag_stdv_infl_fcsts ;
    counts = counts + nb_fcsts ;
    diagM.block(counts,0,nb_fcsts,1)     = diag_stdv_str_fcsts ;
    counts = counts + nb_fcsts ;
    diagM.block(counts,0,nb_fcsts,1)     = diag_stdv_y10_fcsts ;
    M = diagM.asDiagonal() ;
    
    // Predicting the observables
    Obs_pred  = RES_measurement("fX") ;
    M_pred    = G * P_pred * G.transpose() + M * M.transpose();
    
    // if(t == 1){
    //   Rcout << "W_pred : " << W_pred << "\n";
    //   Rcout << "G : " << G << "\n";
    //   Rcout << "Obs_pred : " << Obs_pred << "\n";
    //   Rcout << "Yt : " << Yt << "\n";
    //   Rcout << "X_t : " << X_t << "\n";
    // }
    
    //Rcout << "Obs_pred : " << Obs_pred << "\n";
    
    // UPDATING STEP
    //--------------
    residual = Yt - Obs_pred ;
    
    // Determining which components are observed
    Array<bool,Dynamic,1> bool_na(Yt.array() == Yt.array());
    
    int nb_obs_modif = bool_na.count();
    
    // Fill the new components
    Eigen::MatrixXd resid_modif = Eigen::MatrixXd::Zero(nb_obs_modif, 1);
    Eigen::MatrixXd M_modif     = Eigen::MatrixXd::Zero(nb_obs_modif, nb_obs_modif);
    Eigen::MatrixXd G_modif     = Eigen::MatrixXd::Zero(nb_obs_modif, n_x);
    
    // Control for unobserved components
    //----------------------------------
    int counter_rows = 0;
    for (int i = 0; i < nb_obs; i++){
      
      if(bool_na(i)==TRUE){
        
        resid_modif(counter_rows)  = residual(i);
        G_modif.row(counter_rows)  = G.row(i);
        
        int counter_cols = 0 ;
        
        for (int j = 0; j < nb_obs; j++){
          
          if (bool_na(j) == TRUE){
            
            M_modif(counter_rows, counter_cols) = M_pred(i,j);
            counter_cols += 1;
          }
        }
        counter_rows += 1;
      }
      
    }// End of the NA loop
    
    
    // Perform the update step and the loglik
    //---------------------------------------
    double det_M_pred      = M_modif.determinant();
    
    Eigen::MatrixXd M_inv = M_modif.inverse();
    Gain  = P_pred * G_modif.transpose() * M_inv;
    
    W_up  = W_pred + Gain * resid_modif ;
    P_up  = P_pred - Gain * G_modif * P_pred ;
    
    // if(t == 0){
    //   Rcout << "resid_modif : " << resid_modif << "\n";
    //   Rcout << "M_inv : " << M_inv << "\n";
    //   Rcout << "ZRZ : " << M * M.transpose() << "\n";
    //   Rcout << "GPG : " << G * P_pred * G.transpose() << "\n";
    //   Rcout << "sum((P_pred)) : " << P_pred.sum() << "\n";
    // }

    double lik_val         = -.5*(nb_obs_modif * log(2*pi) + log(det_M_pred) +
                                  (resid_modif.transpose() * M_inv * resid_modif)(0,0));
    
    // Update the values before looping back
    //--------------------------------------
    W = W_up;
    P = P_up;
    
    // if(t == 3){
    //   Rcout << "P0 : " << P0 << "\n";
    //   Rcout << "W_pred : " << W_pred << "\n";
    //   Rcout << "P_pred : " << P_pred << "\n";
    //   Rcout << "W_up : " << W_up << "\n";
    //   Rcout << "P_up : " << P_up << "\n";
    //   Rcout << "Gain : " << Gain << "\n";
    //   Rcout << "residual : " << residual << "\n";
    //   Rcout << "Obs_pred : " << Obs_pred << "\n";
    //   Rcout << "M * M.transpose() : " << M * M.transpose() << "\n";
    // }

    // STORING THE VALUES
    //---------------------
    W_t_t.col(t)          = W_up;
    P_t_t(t)              = P_up;
    
    if(indic_compute_fitted==1){
      RES_measurement = f_measurement(model_sol,
                                      RES_aux,
                                      RES_fcsts,
                                      maturity_nomi_yds,
                                      maturity_real_yds,
                                      horizon_fcsts,
                                      W_up.transpose(),
                                      i_bar) ;
      Obs_fitted_t      = RES_measurement("fX") ;
      Obs_fitted.col(t) = Obs_fitted_t ;
    }
    
    loglik(t) = lik_val;
    
  }
  
  double logl = loglik.sum();
  
  // Sending results back
  //---------------------
  return List::create(
    Named("loglik")         = logl,
    Named("loglik.vector")  = loglik,
    Named("W.updated")      = W_t_t,
    Named("P.updated")      = P_t_t,
    Named("Obs.fitted")     = Obs_fitted,
    Named("list.G")         = list_G,
    Named("M")              = M,
    Named("Y")              = Y);
}


