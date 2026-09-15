// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
using namespace Rcpp;

static double norm_cdf(double x) {
  return R::pnorm5(x, 0.0, 1.0, 1, 0);
}

static double norm_pdf(double x) {
  return R::dnorm4(x, 0.0, 1.0, 0);
}

// [[Rcpp::export]]
Rcpp::List y_fitting_jfec_cpp(const NumericMatrix& X0,
                              const NumericVector& A_X_for,
                              const NumericMatrix& B_X_for,
                              const NumericVector& A_X_exp,
                              const double r_lb,
                              const NumericVector& s_n,
                              const NumericVector& A_X_for_pi,
                              const NumericMatrix& B_X_pi,
                              const NumericMatrix& Sigma2_X,
                              const NumericMatrix& B_X_cum,
                              const NumericMatrix& B_X_cum_pi) {
  int K = X0.nrow();
  int T0 = X0.ncol();
  int Mm = A_X_for.size();

  NumericMatrix yfit_all_n(T0, Mm);
  NumericMatrix yfit_all_r(T0, Mm);
  NumericMatrix JJ_n(Mm, K * T0);
  NumericMatrix JJ_r(Mm, K * T0);
  NumericMatrix Probs(T0, Mm);

  NumericVector BXcumSig2BXcum_pi(Mm);
  BXcumSig2BXcum_pi[0] = 0.0;
  for (int m = 1; m < Mm; ++m) {
    double val = 0.0;
    for (int i = 0; i < K; ++i) {
      double sig_bcum_pi = 0.0;
      for (int j = 0; j < K; ++j) {
        sig_bcum_pi += Sigma2_X(i, j) * B_X_cum_pi(j, m);
      }
      val += B_X_cum(i, m) * sig_bcum_pi;
    }
    BXcumSig2BXcum_pi[m] = val;
  }

  for (int t = 0; t < T0; ++t) {
    NumericVector f_fit(Mm);
    NumericVector r_fit(Mm);

    double shadow_short = A_X_exp[0];
    for (int k = 0; k < K; ++k) {
      shadow_short += X0(k, t) * B_X_for(k, 0);
    }
    double prob_short = 1.0;
    double y_short = shadow_short;
    if (y_short < r_lb) {
      prob_short = 0.0;
      y_short = r_lb;
    }
    Probs(t, 0) = prob_short;
    f_fit[0] = y_short;

    double pi_fit_0 = A_X_for_pi[0];
    for (int k = 0; k < K; ++k) {
      pi_fit_0 += X0(k, t) * B_X_pi(k, 0);
    }
    r_fit[0] = f_fit[0] - pi_fit_0;

    for (int k = 0; k < K; ++k) {
      int row = t * K + k;
      double jj_f = (shadow_short > r_lb) ? B_X_for(k, 0) : 0.0;
      JJ_n(0, row) = jj_f;
      JJ_r(0, row) = jj_f - B_X_pi(k, 0);
    }

    for (int m = 1; m < Mm; ++m) {
      double mu = A_X_exp[m] - r_lb;
      for (int k = 0; k < K; ++k) {
        mu += X0(k, t) * B_X_for(k, m);
      }
      double z_n = mu / s_n[m - 1];
      double prob_long = norm_cdf(z_n);
      double pdf_z = norm_pdf(z_n);
      Probs(t, m) = prob_long;

      // Evaluate the unexpanded Wu-Xia term, with convexity inside g.
      // The real-rate interaction below uses the unshifted probability.
      double nominal_mu = mu + A_X_for[m] - A_X_exp[m];
      double nominal_z = nominal_mu / s_n[m - 1];
      double nominal_probability = norm_cdf(nominal_z);
      f_fit[m] = r_lb + nominal_mu * nominal_probability +
        norm_pdf(nominal_z) * s_n[m - 1];

      double pi_fit = A_X_for_pi[m];
      for (int k = 0; k < K; ++k) {
        pi_fit += X0(k, t) * B_X_pi(k, m);
      }
      r_fit[m] = f_fit[m] - pi_fit + prob_long * BXcumSig2BXcum_pi[m];

      for (int k = 0; k < K; ++k) {
        int row = t * K + k;
        double jj_f = nominal_probability * B_X_for(k, m);
        double jj_real = jj_f - B_X_pi(k, m) +
          (pdf_z / s_n[m - 1]) * B_X_for(k, m) * BXcumSig2BXcum_pi[m];
        JJ_n(m, row) = JJ_n(m - 1, row) + jj_f;
        JJ_r(m, row) = JJ_r(m - 1, row) + jj_real;
      }
    }

    double cum_n = 0.0;
    double cum_r = 0.0;
    for (int m = 0; m < Mm; ++m) {
      cum_n += f_fit[m];
      cum_r += r_fit[m];
      yfit_all_n(t, m) = cum_n / (m + 1);
      yfit_all_r(t, m) = cum_r / (m + 1);
      for (int k = 0; k < K; ++k) {
        int row = t * K + k;
        JJ_n(m, row) /= (m + 1);
        JJ_r(m, row) /= (m + 1);
      }
    }
  }

  return List::create(
    _["yfit_all_n"] = yfit_all_n,
    _["yfit_all_r"] = yfit_all_r,
    _["JJ_n"] = JJ_n,
    _["JJ_r"] = JJ_r,
    _["Probs"] = Probs
  );
}
