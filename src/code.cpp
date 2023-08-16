#include <Rcpp.h>
#include <RcppEigen.h>
#include "TVopt_custom.h"

// #include <RcppParallel.h>

using namespace Rcpp;
// using namespace RcppParallel;

// Linear Algebra helper code
typedef Eigen::SparseMatrix<double> SpMat;

SpMat difference_matrix(int n) {
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(2 * (n - 1));

  for (int i = 0; i < n - 1; i++) {
    tripletList.push_back(Eigen::Triplet<double>(i, i, -1.0));
    tripletList.push_back(Eigen::Triplet<double>(i, i + 1, 1.0));
  }

  SpMat diff(n - 1, n);
  diff.setFromTriplets(tripletList.begin(), tripletList.end());

  return diff;
}



// Total variation penalization
//' @export
// [[Rcpp::export]]
NumericVector Rcpp_classicTautString_TV1(NumericVector input, double lambda) {
  int length = input.size();
  NumericVector output(length);
  classicTautString_TV1(input.begin(), length, lambda, output.begin());
  return output;
}


int searchsortedfirst(const std::vector<double>& arr, double target) {
  for (size_t i = 0; i < arr.size(); i++) {
    if (arr[i] <= target) {
      return i;
    }
  }
  return arr.size(); // Return size of the array if not found
}

class Grenander {
public:
  NumericVector x_knots;
  NumericVector y_knots;
  NumericVector slope_knots;
  int length;
  std::vector<double> shifted_fs_left;
  std::vector<double> shifted_fs_right;

  Grenander(NumericVector x_knots, NumericVector y_knots, NumericVector slope_knots) :
    x_knots(x_knots), y_knots(y_knots), slope_knots(slope_knots),
    length(x_knots.size()),
    shifted_fs_left(slope_knots.size()), shifted_fs_right(slope_knots.size())
  {
    // Check if the sizes are consistent
    if (y_knots.size() != x_knots.size() || slope_knots.size() != (x_knots.size() - 1)) {
      stop("Inconsistent size of input vectors");
    }
  }

  // Calculate the CDF at point t
  double cdf(double t) {
    if (t >= x_knots[length - 1]) {
      return y_knots[length - 1];
    } else if (t <= x_knots[0]) {
      return 0.0;
    } else {
      int idx = std::upper_bound(x_knots.begin(), x_knots.end(), t) - x_knots.begin() - 1;

      if (idx < 0 || idx >= length - 1) { // Safeguard against out-of-bounds access
        stop("Invalid index detected");
      }

      double loc_L = x_knots[idx];
      double F_L = y_knots[idx];
      double loc_R = x_knots[idx+1];
      double F_R = y_knots[idx+1];
      double lambda = (t - loc_L) / (loc_R - loc_L);
      return (1 - lambda) * F_L + lambda * F_R;
    }
  }

  double pdf(double t) {
    if(t > x_knots[length-1]) {
      return 0.0;
    } else if (t <= x_knots[0]) {  // Change condition here to include the equality
      return slope_knots[0];
    } else {
      int idx = std::lower_bound(x_knots.begin(), x_knots.end(), t) - x_knots.begin();
      idx -= 1;
      return slope_knots[idx];
    }
  }

  std::pair<double, double> invert_subgradient(double lambda) {
    double loc_lb, loc_ub;

    // λ is greater than the first value of fs (slope_knots in C++)
    if (lambda > slope_knots[0]) {
      loc_lb = 0.0;
      loc_ub = 0.0;
    }
    // λ is less than the last value of fs (slope_knots in C++)
    else if (lambda < slope_knots[length - 2]) { // adjusted this line
      loc_lb = 1.0;
      loc_ub = 1.0;
    }
    // Find the value in fs (slope_knots in C++) just greater than or equal to λ
    else {
      auto comp = [](double a, double b) { return a > b; };
      int idx = std::lower_bound(slope_knots.begin(), slope_knots.end(), lambda, comp) - slope_knots.begin();

      loc_lb = x_knots[idx];
      double loc_lb_p1 = x_knots[idx + 1];
      loc_ub = (this->pdf(loc_lb_p1) == lambda) ? loc_lb_p1 : loc_lb;
    }

    return std::make_pair(loc_lb, loc_ub);
  }

  double invert_subgradient_tilted(double mu, double rho, double b, double alpha) {
    double mu_n = mu / (1 + alpha * mu);
    double rho_n = rho / (1 + alpha * mu);
    double mu_p = mu_n - rho_n * b;

    // Assume shifted_fs_left and shifted_fs_right are declared as member variables or elsewhere
    for (int i = 0; i < length - 1; i++) {
      shifted_fs_left[i] = slope_knots[i] - rho_n * x_knots[i];
      shifted_fs_right[i] = slope_knots[i] - rho_n * x_knots[i + 1];
    }

    double loc;
    if (mu_p >= shifted_fs_left[0]) {
      loc = 0.0;
    } else if (mu_p < shifted_fs_right[length - 2]) {  // adjusted index for last element
      loc = 1.0;
    } else {
      auto comp = [](double a, double b) { return a > b; };
      auto idx = std::upper_bound(shifted_fs_left.begin(), shifted_fs_left.end(), mu_p, comp);
      int index = std::distance(shifted_fs_left.begin(), idx);

      if (index >= 1 && shifted_fs_left[index-1] <= mu_p) {
        index = index - 1;
      }

      double loc_lb = x_knots[index];
      double loc_lb_m1 = (index > 0) ? x_knots[index - 1] : 0.0;  // guarding against negative indexing
      if (shifted_fs_right[index - 1] >= mu_p) {
        loc = loc_lb;
      } else {
        double alpha_ = (mu_p - slope_knots[index - 1] + rho_n * loc_lb) / (rho_n * (loc_lb - loc_lb_m1));
        loc = alpha_ * loc_lb_m1 + (1 - alpha_) * loc_lb;
      }
    }
    return loc;
  }


};

//' @export
// [[Rcpp::export]]
double Grenander_cdf_export(List grenander_list, double t) {
  // Convert list to Grenander
  NumericVector x_knots = grenander_list["locs"];
  NumericVector y_knots = grenander_list["Fs"];
  NumericVector slope_knots = grenander_list["fs"];

  Grenander grenander(x_knots, y_knots, slope_knots);

  // Compute the CDF
  return grenander.cdf(t);
}

//' @export
// [[Rcpp::export]]
double Grenander_pdf_export(List grenander_list, double t) {
  // Convert list to Grenander
  NumericVector x_knots = grenander_list["locs"];
  NumericVector y_knots = grenander_list["Fs"];
  NumericVector slope_knots = grenander_list["fs"];

  Grenander grenander(x_knots, y_knots, slope_knots);

  // Compute the PDF
  return grenander.pdf(t);
}

//' @export
// [[Rcpp::export]]
NumericVector Grenander_invert_subgradient(List grenander_list, double lambda) {
  NumericVector x_knots = grenander_list["locs"];
  NumericVector y_knots = grenander_list["Fs"];
  NumericVector slope_knots = grenander_list["fs"];

  Grenander grenander(x_knots, y_knots, slope_knots);

  std::pair<double, double> interval = grenander.invert_subgradient(lambda);

  return NumericVector::create(interval.first, interval.second);
}


template <typename Func>
double bisection(double a, double b, double tol, Func&& fun, bool use_geometric_mean = false) {
  double t;
  int max_iterations = 100000;  // or another appropriate value
  int iterations = 0;
  double fb = fun(b);
  double fa = fun(a);
  double ft;

  while (fb - fa > tol) {

    if (iterations++ > max_iterations) {
      std::cout << "Warning: Maximum number of iterations reached!" << std::endl;
      break;
    }

    if (use_geometric_mean) {
      t = std::sqrt(a * b);  // Geometric mean for positive a and b
    } else {
      t = (a + b) / 2;  // Arithmetic mean
    }
    ft = fun(t);
    if (ft < 0) {
      a = t;
      fa = ft;
    } else {
      b = t;
      fb = ft;
    }
  }
  return (a + b) / 2;
}



class GrenMix {
public:
  std::vector<Grenander*> grenanders; // Store pointers to Grenander objects
  NumericVector test_ms;
  std::vector<double> ts_values;
  std::vector<double> Fs_values;

  GrenMix(List train_grenander, NumericVector test_ms_input) {
    NumericVector test_ms_copy = Rcpp::clone(test_ms_input);
    double test_ms_sum = sum(test_ms_copy);
    test_ms_copy = test_ms_copy / test_ms_sum; // Normalize test_ms_copy

    this->test_ms = test_ms_copy;


    for (int i = 0; i < train_grenander.size(); i++) {
      List grenander = train_grenander[i];
      Grenander *g = new Grenander(
        as<NumericVector>(grenander["locs"]),
        as<NumericVector>(grenander["Fs"]),
        as<NumericVector>(grenander["fs"])
      );
      grenanders.push_back(g);
    }

    ts_values.resize(grenanders.size());
    Fs_values.resize(grenanders.size());
  }

  double pdf(double t) {
    double total_pdf = 0.0;
    for (int i = 0; i < grenanders.size(); i++) {
      total_pdf += test_ms[i] * grenanders[i]->pdf(t);
    }
    return total_pdf;
  }

  double cdf(double t) {
    double total_cdf = 0.0;
    for (int i = 0; i < grenanders.size(); i++) {
      total_cdf += test_ms[i] * grenanders[i]->cdf(t);
    }
    return total_cdf;
  }

  std::pair<double, double> lagrange_balance(double lambda, double alpha, bool linearized) {
    double sum_ts_L = 0.0, sum_ts_R = 0.0, sum_Fs_L = 0.0, sum_Fs_R = 0.0;

    for (int i = 0; i < grenanders.size(); i++) {
      std::pair<double, double> ts = grenanders[i]->invert_subgradient(lambda);
      double ts_L = ts.first;
      double ts_R = ts.second;

      double Fs_L = grenanders[i]->cdf(ts_L);
      double Fs_R = grenanders[i]->cdf(ts_R);

      sum_ts_L += ts_L * test_ms[i];
      sum_ts_R += ts_R * test_ms[i];
      sum_Fs_L += Fs_L * test_ms[i];
      sum_Fs_R += Fs_R * test_ms[i];
    }

    double Fdr_L, Fdr_R;
    if (linearized) {
      Fdr_L = (sum_ts_L - alpha * sum_Fs_L);
      Fdr_R = (sum_ts_R - alpha * sum_Fs_R);
    } else {
      Fdr_L = (sum_Fs_L == 0.0 ? 0.0 : sum_ts_L / sum_Fs_L);
      Fdr_R = (sum_Fs_R == 0.0 ? 0.0 : sum_ts_R / sum_Fs_R);
    }

    return std::make_pair(Fdr_L, Fdr_R);
  }

  double lagrange_balance_tilted(double mu, const std::vector<double>& bs, const std::vector<double>& rhos, double alpha, bool linearized) {
    double sum_ts = 0.0, sum_Fs = 0.0;

    for (int i = 0; i < grenanders.size(); i++) {
      ts_values[i] = grenanders[i]->invert_subgradient_tilted(mu, rhos[i], bs[i], alpha);
      Fs_values[i] = grenanders[i]->cdf(ts_values[i]);

      sum_ts += ts_values[i] * test_ms[i];
      sum_Fs += Fs_values[i] * test_ms[i];
    }

    double Fdr;
    if (linearized) {
      Fdr = (sum_ts - alpha * sum_Fs);
    } else {
      Fdr = (sum_Fs == 0.0 ? 0.0 : sum_ts / sum_Fs);
    }

    return Fdr;
  }

  std::vector<double> unregularized_thresholds_bh(double alpha) {
    std::vector<double> all_fs;
    for (size_t i = 0; i < grenanders.size(); i++) {
      std::vector<double> fs(grenanders[i]->slope_knots.begin(), grenanders[i]->slope_knots.end());
      all_fs.insert(all_fs.end(), fs.begin(), fs.end());
    }

    std::sort(all_fs.rbegin(), all_fs.rend()); // sort in descending order

    double lambda_opt;
    for (size_t i = 0; i < all_fs.size(); i++) {
      double lambda = all_fs[i];
      std::pair<double, double> balance = lagrange_balance(lambda, alpha, false);
      if (alpha >= balance.first && alpha <= balance.second) {
        lambda_opt = lambda;
        break;
      }
    }

    std::vector<double> ts_L, ts_R;
    for (size_t i = 0; i < grenanders.size(); i++) {
      std::pair<double, double> ts = grenanders[i]->invert_subgradient(lambda_opt);
      ts_L.push_back(ts.first);
      ts_R.push_back(ts.second);
    }


    std::pair<double, double> lin_Fdr = lagrange_balance(lambda_opt, alpha, true);
    double tmp_comb = lin_Fdr.second / (lin_Fdr.second - lin_Fdr.first);

    std::vector<double> result;
    for (size_t i = 0; i < ts_L.size(); i++) {
      result.push_back(tmp_comb * ts_L[i] + (1 - tmp_comb) * ts_R[i]);
    }

    return result;
  }

  std::pair<double, double> single_t_FDR(double alpha); // defined below
  double lambda_max(double alpha, double mu_dual, double single_t); // defined below
  NumericMatrix total_variation_path_ts(double alpha, const NumericVector& lambda_multipliers); // defined below

  ~GrenMix() {
    for (std::vector<Grenander*>::iterator it = grenanders.begin(); it != grenanders.end(); ++it) {
      delete *it;
    }
  }
};


struct FdrFunctor {
  GrenMix& mix;
  double alpha;

  FdrFunctor(GrenMix& mix, double alpha) : mix(mix), alpha(alpha) {}

  double operator()(double t) {
    return t == 0 ? 1 / mix.pdf(t) - alpha : t / mix.cdf(t) - alpha;
  }
};


std::pair<double, double> GrenMix::single_t_FDR(double alpha) {
  double t, mu_dual;
  if (pdf(0.0) <= 1 / alpha) {
    t = 0.0;
  } else {
    double a = 0.0, b = 1.0;
    double tol = alpha * std::pow(10.0, -4.0);  // Set tolerance.
    FdrFunctor fun(*this, alpha);
    t = bisection(a, b, tol, fun);
  }
  mu_dual = pdf(t) / (1 - alpha * pdf(t));
  return std::make_pair(t, mu_dual);
}

// [[Rcpp::depends(RcppEigen)]]
double GrenMix::lambda_max(double alpha, double mu_dual, double single_t) {
  std::vector<double> pdf_values;
  for (std::vector<Grenander*>::iterator it = grenanders.begin(); it != grenanders.end(); ++it) {
    pdf_values.push_back((*it)->pdf(single_t));
  }


  Eigen::Map<Eigen::VectorXd> pdf_map(pdf_values.data(), pdf_values.size());
  Eigen::VectorXd gren_sub = pdf_map.array() * (1 + alpha * mu_dual) - mu_dual;

  Eigen::Map<Eigen::VectorXd> test_ms_map = as<Eigen::Map<Eigen::VectorXd>>(test_ms);
  Eigen::VectorXd ys = test_ms_map.array() * gren_sub.array();

  SpMat D = difference_matrix(ys.size());

  SpMat DD = D * D.transpose();

  // Define the solver
  Eigen::SparseLU<SpMat> solver;
  solver.compute(DD);

  if(solver.info()!=Eigen::Success) {
    // Decomposition failed
    throw std::runtime_error("Decomposition Failed");
  }

  Eigen::VectorXd lambda = solver.solve(D * ys);

  if(solver.info()!=Eigen::Success) {
    // Solving failed
    throw std::runtime_error("Solving Failed");
  }

  double lambda_max = lambda.cwiseAbs().maxCoeff();

  return lambda_max;
}

NumericMatrix GrenMix::total_variation_path_ts(double alpha, const NumericVector& lambda_multipliers) {
  int num_grenanders = this->grenanders.size();
  NumericMatrix ts_array(num_grenanders, lambda_multipliers.size()+2);

  double max_density = 0.0;

  for(int i = 0; i < num_grenanders; ++i) {
    double first_density = this->grenanders[i]->slope_knots[0];
    if(first_density > max_density) {
      max_density = first_density;
    }
  }

  // std::cout << "max_density: " << max_density << std::endl;

  if (max_density < 1/alpha) {
    return ts_array;
  }

  double min_density = max_density;
  for(int i = 0; i < num_grenanders; ++i) {
    double last_density = this->grenanders[i]->slope_knots[this->grenanders[i]->slope_knots.size() - 1];  // Get the last element
    if(last_density < min_density && last_density > 0.0) {
      min_density = last_density;
    }
  }
  // std::cout << "Just before single T " << std::endl;

  std::pair<double, double> mu_t = this->single_t_FDR(alpha);
  double single_t = mu_t.first;
  double mu_dual = mu_t.second;

  // std::cout << "Just after single T " << std::endl;


  double lambda_max = this->lambda_max(alpha, mu_dual, single_t);

  //std::cout << "lambda_max " << lambda_max << std::endl;

  // Create lambdas using lambda_multipliers
  int num_lambdas = lambda_multipliers.size();
  std::vector<double> lambdas(num_lambdas);
  for(int i = 0; i < num_lambdas; ++i) {
    lambdas[i] = lambda_max * lambda_multipliers[i];
  }

  std::vector<double> bs(num_grenanders, single_t);
  std::vector<double> us(num_grenanders, 0.0);
  std::vector<double> ts_unreg = this->unregularized_thresholds_bh(alpha);
  //std::vector<double> ts_unreg(num_grenanders, 0.1);


  for(int j = 0; j < num_grenanders; ++j){
    ts_array(j, 0) = single_t;
    ts_array(j, lambdas.size()+1) = ts_unreg[j];
  }

  // preallocate outside the loop
  std::vector<double> rho_prop(num_grenanders, 0.0);
  std::vector<double> ts_iter(num_grenanders, 0.0);
  std::vector<double> shift_iter_temp(num_grenanders, 0.0);
  std::vector<double> previous_shift_iter_temp(num_grenanders, 10.0);

  double tol = alpha * std::pow(10.0, -4.0);  // Set tolerance.

  for(int i = 0; i < lambdas.size(); ++i) {
    double rho = 100 * lambdas[i];// * 5000;
    double lambda_by_rho = lambdas[i] / rho;
    //if (i > 0){
    //  double mult = lambda_multipliers[i-1]/lambda_multipliers[i];

    //  for (int k = 0; k < num_grenanders; k++) {
    //    us[k] = us[k] * mult;
    //  }
    //}

    for (int k = 0; k < num_grenanders; k++) {
      rho_prop[k] = rho / this->test_ms[k];
    }

    //std::cout << "lambda: " << rho << std::endl;

    for (int k = 0; k < num_grenanders; k++) {
      rho_prop[k] = rho / this->test_ms[k];
    }

    for(int j = 0; j < 500; ++j) { // should replace this with an actual stopping rule

      if (i > 0){
        // break;
      }

      double zero_val = bisection(min_density, max_density, tol, [this, &bs, &rho_prop, &alpha](double mu) -> double {
        return -(this->lagrange_balance_tilted(mu, bs, rho_prop, alpha, false) - alpha);
      }, true);

      //std::cout << "zero_val " << zero_val << std::endl;
      //std::cout << "FDR " << lagrange_balance_tilted(zero_val, bs, rho_prop, alpha, false) << std::endl;

      for (int k = 0; k < num_grenanders; k++) {
        ts_iter[k] = this->grenanders[k]->invert_subgradient_tilted(zero_val, rho_prop[k], bs[k], alpha);
      }

      //std::cout << "ts_iter0 " << ts_iter[0] << std::endl;


      for (int k = 0; k < num_grenanders; k++) {
        bs[k] =  ts_iter[k] + us[k];
      }

      classicTautString_TV1(bs.data(), num_grenanders, lambda_by_rho, shift_iter_temp.data());
      for (double &value : shift_iter_temp) {
        if (value < 0.0) {
          value = 0.0;
        } else if (value > 1.0) {
          value = 1.0;
        }
      }

      //std::cout << "shift_iter " << shift_iter_temp[0] << std::endl;

      double diff_sum = 0.0;
      double prev_sum = 0.0;
      for (int k = 0; k < num_grenanders; k++) {
        diff_sum += std::abs(shift_iter_temp[k] - ts_iter[k]);
        prev_sum += std::abs(ts_iter[k]);
      }

      double relative_change = (prev_sum > 0) ? (diff_sum / prev_sum) : 0;

      if (relative_change < 1e-3) {
         break;  // break out of the loop if relative change is below the threshold
      }


      for (int k = 0; k < num_grenanders; k++) {
        us[k] = ts_iter[k] + us[k] - shift_iter_temp[k];
        bs[k] = shift_iter_temp[k] - us[k];
        previous_shift_iter_temp[k] = shift_iter_temp[k];
      }
      //std::cout << "us0 " << us[0] << std::endl;

    }
    for(int j = 0; j < num_grenanders; ++j){
      ts_array(j,i+1) = shift_iter_temp[j];
    }
  }
  return ts_array;
}



 // [[Rcpp::export]]
 NumericVector single_t_FDR(List train_grenander, NumericVector test_ms, double alpha) {
   GrenMix gren_mix(train_grenander, test_ms);
   std::pair<double, double> result = gren_mix.single_t_FDR(alpha);
   return NumericVector::create(Named("t") = result.first, Named("mu_dual") = result.second);
 }

// [[Rcpp::export]]
double GrenMix_pdf(List train_grenander, NumericVector test_ms, double t) {
  GrenMix gren_mix(train_grenander, test_ms);
  return gren_mix.pdf(t);
}

// [[Rcpp::export]]
double GrenMix_cdf(List train_grenander, NumericVector test_ms, double t) {
  GrenMix gren_mix(train_grenander, test_ms);
  return gren_mix.cdf(t);
}

// [[Rcpp::export]]
Rcpp::NumericVector GrenMix_lagrange_balance(Rcpp::List train_grenander, Rcpp::NumericVector test_ms,
                                             double lambda, double alpha, bool linearized) {
  GrenMix grenMix(train_grenander, test_ms);
  std::pair<double, double> result = grenMix.lagrange_balance(lambda, alpha, linearized);
  return Rcpp::NumericVector::create(result.first, result.second);
}

 // [[Rcpp::export]]
 NumericVector unregularized_thresholds_bh(List train_grenander, NumericVector test_ms, double alpha) {
   GrenMix grenmix(train_grenander, test_ms);
   std::vector<double> result = grenmix.unregularized_thresholds_bh(alpha);
   return Rcpp::wrap(result);
 }

// [[Rcpp::export]]
double InvertSubgradientTilted(List grenander, double mu, double rho, double b, double alpha) {
  Grenander *g = new Grenander(
    as<NumericVector>(grenander["locs"]),
    as<NumericVector>(grenander["Fs"]),
    as<NumericVector>(grenander["fs"])
  );

  return g->invert_subgradient_tilted(mu, rho, b, alpha);
}

// [[Rcpp::export]]
double lambdaMax(List train_grenander, NumericVector test_ms, double alpha, double mu_dual, double single_t) {
  GrenMix gren_mix(train_grenander, test_ms);
  return gren_mix.lambda_max(alpha, mu_dual, single_t);
}

 // [[Rcpp::export]]
 NumericMatrix optimal_ts(List train_grenander, NumericVector test_ms, double alpha, NumericVector lambda_multipliers) {
   GrenMix gren_mix(train_grenander, test_ms);
   return gren_mix.total_variation_path_ts(alpha, lambda_multipliers);
 }


