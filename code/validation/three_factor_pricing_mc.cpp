// Independent Monte Carlo benchmark for the three-factor approximation study.
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List three_factor_bond_mc(const NumericMatrix& states,
                                const NumericVector& rho,
                                const NumericVector& innovation_sd,
                                const NumericVector& shadow_loading,
                                const NumericVector& inflation_loading,
                                const double shadow_intercept,
                                const double inflation_intercept,
                                const double lower_bound,
                                const IntegerVector& horizons,
                                const NumericMatrix& control_price_nominal,
                                const NumericMatrix& control_price_real,
                                const int antithetic_pairs,
                                const int seed) {
  const int ns = states.nrow();
  const int nh = horizons.size();
  const int hmax = max(horizons);
  const int npaths = 2 * antithetic_pairs;
  NumericMatrix sum_dn(ns, nh), sumsq_dn(ns, nh), sum_dr(ns, nh), sumsq_dr(ns, nh);
  std::vector<int> hmap(hmax + 1, -1);
  for (int j = 0; j < nh; ++j) hmap[horizons[j]] = j;

  Function set_seed("set.seed");
  set_seed(seed);
  RNGScope scope;

  // Process one antithetic pair at a time. This deliberately does not use any
  // analytical-approximation objects or routines.
  for (int draw = 0; draw < antithetic_pairs; ++draw) {
    NumericMatrix shock(hmax, 3);
    NumericMatrix pair_dn(ns, nh), pair_dr(ns, nh);
    for (int h = 0; h < hmax; ++h)
      for (int k = 0; k < 3; ++k)
        shock(h, k) = R::rnorm(0.0, innovation_sd[k]);

    for (int sign_index = 0; sign_index < 2; ++sign_index) {
      const double sg = sign_index == 0 ? 1.0 : -1.0;
      for (int s = 0; s < ns; ++s) {
        double x[3] = {states(s, 0), states(s, 1), states(s, 2)};
        double log_n = 0.0, log_r = 0.0, log_un_n = 0.0, log_un_r = 0.0;
        for (int h = 1; h <= hmax; ++h) {
          double shadow = shadow_intercept;
          for (int k = 0; k < 3; ++k) shadow += shadow_loading[k] * x[k];
          const double effective_rate = std::max(lower_bound, shadow);
          log_n -= effective_rate;
          log_r -= effective_rate;
          log_un_n -= shadow;
          log_un_r -= shadow;

          for (int k = 0; k < 3; ++k)
            x[k] = rho[k] * x[k] + sg * shock(h - 1, k);
          double inflation = inflation_intercept;
          for (int k = 0; k < 3; ++k) inflation += inflation_loading[k] * x[k];
          log_r += inflation;
          log_un_r += inflation;

          const int col = hmap[h];
          if (col >= 0) {
            const double pn = std::exp(log_n);
            const double pr = std::exp(log_r);
            const double pun = std::exp(log_un_n);
            const double pur = std::exp(log_un_r);
            pair_dn(s, col) += 0.5 * (pn - pun);
            pair_dr(s, col) += 0.5 * (pr - pur);
          }
        }
      }
    }
    for (int s = 0; s < ns; ++s) for (int j = 0; j < nh; ++j) {
      sum_dn(s, j) += pair_dn(s, j); sumsq_dn(s, j) += pair_dn(s, j) * pair_dn(s, j);
      sum_dr(s, j) += pair_dr(s, j); sumsq_dr(s, j) += pair_dr(s, j) * pair_dr(s, j);
    }
  }

  NumericMatrix price_n(ns, nh), price_r(ns, nh), se_n(ns, nh), se_r(ns, nh);
  for (int s = 0; s < ns; ++s) for (int j = 0; j < nh; ++j) {
    const double mdn = sum_dn(s, j) / antithetic_pairs;
    const double mdr = sum_dr(s, j) / antithetic_pairs;
    price_n(s, j) = control_price_nominal(s, j) + mdn;
    price_r(s, j) = control_price_real(s, j) + mdr;
    const double vn = std::max(0.0, sumsq_dn(s, j) / antithetic_pairs - mdn * mdn);
    const double vr = std::max(0.0, sumsq_dr(s, j) / antithetic_pairs - mdr * mdr);
    se_n(s, j) = std::sqrt(vn / antithetic_pairs);
    se_r(s, j) = std::sqrt(vr / antithetic_pairs);
  }
  return List::create(_["price_nominal"] = price_n,
                      _["price_real"] = price_r,
                      _["se_price_nominal"] = se_n,
                      _["se_price_real"] = se_r,
                      _["npaths"] = npaths);
}
