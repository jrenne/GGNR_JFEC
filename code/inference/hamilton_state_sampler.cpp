// Draw pointwise filtered states from the Gaussian approximation delivered by
// the reduced-form EKF. Random normals are generated in R so results remain
// reproducible under parallel execution.
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

typedef std::vector<double> DenseMatrix;

inline double& at(DenseMatrix& value, int row, int column, int ncol) {
  return value[row * ncol + column];
}
inline double atc(const DenseMatrix& value, int row, int column, int ncol) {
  return value[row * ncol + column];
}

DenseMatrix covariance_slice(const NumericVector& values, int time,
                             int dimension) {
  DenseMatrix output(dimension * dimension);
  int offset = dimension * dimension * time;
  for (int column = 0; column < dimension; ++column) {
    for (int row = 0; row < dimension; ++row) {
      at(output, row, column, dimension) =
        values[offset + row + dimension * column];
    }
  }
  return output;
}

DenseMatrix cholesky_with_jitter(const DenseMatrix& covariance,
                                 int dimension) {
  double scale = 1.0;
  for (int index = 0; index < dimension; ++index) {
    scale = std::max(
      scale, std::abs(atc(covariance, index, index, dimension))
    );
  }
  for (int trial = 0; trial < 8; ++trial) {
    double jitter = scale * 1e-14 * std::pow(10.0, trial);
    DenseMatrix lower(dimension * dimension, 0.0);
    bool valid = true;
    for (int row = 0; row < dimension && valid; ++row) {
      for (int column = 0; column <= row; ++column) {
        double sum = 0.5 * (
          atc(covariance, row, column, dimension) +
            atc(covariance, column, row, dimension)
        );
        if (row == column) sum += jitter;
        for (int inner = 0; inner < column; ++inner) {
          sum -= atc(lower, row, inner, dimension) *
            atc(lower, column, inner, dimension);
        }
        if (row == column) {
          if (!(sum > 0.0) || !std::isfinite(sum)) valid = false;
          else at(lower, row, column, dimension) = std::sqrt(sum);
        } else {
          at(lower, row, column, dimension) = sum /
            atc(lower, column, column, dimension);
        }
      }
    }
    if (valid) return lower;
  }
  stop("Unable to stabilize a filtered-state covariance matrix.");
}

// [[Rcpp::export]]
NumericMatrix draw_filtered_states_cpp(
    const NumericMatrix& state_mean,
    const NumericVector& state_covariance,
    const NumericMatrix& standard_normals) {
  int observation_count = state_mean.nrow();
  int state_count = state_mean.ncol();
  IntegerVector dimensions = state_covariance.attr("dim");
  if (dimensions.size() != 3 || dimensions[0] != state_count ||
      dimensions[1] != state_count ||
      dimensions[2] != observation_count ||
      standard_normals.nrow() != observation_count ||
      standard_normals.ncol() != state_count) {
    stop("Incompatible dimensions in filtered-state sampler.");
  }

  NumericMatrix draws(observation_count, state_count);
  for (int time = 0; time < observation_count; ++time) {
    DenseMatrix lower = cholesky_with_jitter(
      covariance_slice(state_covariance, time, state_count), state_count
    );
    for (int row = 0; row < state_count; ++row) {
      double value = state_mean(time, row);
      for (int inner = 0; inner <= row; ++inner) {
        value += atc(lower, row, inner, state_count) *
          standard_normals(time, inner);
      }
      draws(time, row) = value;
    }
  }
  colnames(draws) = colnames(state_mean);
  return draws;
}
