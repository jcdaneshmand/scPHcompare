#include <Rcpp.h>
#include <algorithm>
#include <functional>
#include <vector>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::List top_k_landscape_profile_cpp(
    Rcpp::NumericMatrix pd, int dimension, Rcpp::NumericVector grid,
    int max_k) {
  if (pd.ncol() < 3 || max_k < 1) {
    Rcpp::stop("pd must have three columns and max_k must be positive.");
  }
  std::vector<double> births;
  std::vector<double> deaths;
  births.reserve(pd.nrow());
  deaths.reserve(pd.nrow());
  for (int i = 0; i < pd.nrow(); ++i) {
    const double dim = pd(i, 0);
    const double birth = pd(i, 1);
    const double death = pd(i, 2);
    if (R_finite(dim) && static_cast<int>(dim) == dimension &&
        R_finite(birth) && R_finite(death) && birth < death) {
      births.push_back(birth);
      deaths.push_back(death);
    }
  }

  Rcpp::NumericMatrix values(grid.size(), max_k);
  Rcpp::NumericVector total_squared(grid.size());
  Rcpp::IntegerVector active_levels(grid.size());
  std::vector<double> tents;
  tents.reserve(births.size());
  for (int g = 0; g < grid.size(); ++g) {
    tents.clear();
    double total = 0.0;
    const double t = grid[g];
    for (std::size_t i = 0; i < births.size(); ++i) {
      const double height = std::min(t - births[i], deaths[i] - t);
      if (height > 0.0) {
        tents.push_back(height);
        total += height * height;
      }
    }
    active_levels[g] = static_cast<int>(tents.size());
    total_squared[g] = total;
    if (tents.size() > static_cast<std::size_t>(max_k)) {
      std::nth_element(
        tents.begin(), tents.begin() + max_k, tents.end(),
        std::greater<double>()
      );
      tents.resize(max_k);
    }
    std::sort(tents.begin(), tents.end(), std::greater<double>());
    for (std::size_t k = 0; k < tents.size(); ++k) values(g, k) = tents[k];
  }
  return Rcpp::List::create(
    Rcpp::Named("values") = values,
    Rcpp::Named("total_squared") = total_squared,
    Rcpp::Named("active_levels") = active_levels,
    Rcpp::Named("finite_intervals") = static_cast<int>(births.size())
  );
}
