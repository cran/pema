// Generated by rstantools.  Do not edit by hand.

/*
    pema is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    pema is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with pema.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by %%NAME%% %%VERSION%%
#include <stan/model/model_header.hpp>
namespace model_prior_meff_namespace {
inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    std::stringstream msg;
    msg << "Found dimension size less than one in simplex declaration"
        << "; variable=" << var_name << "; dimension size expression=" << expr
        << "; expression value=" << val;
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    std::stringstream msg;
    if (val == 1) {
      msg << "Found dimension size one in unit vector declaration."
          << " One-dimensional unit vector is discrete"
          << " but the target distribution must be continuous."
          << " variable=" << var_name << "; dimension size expression=" << expr;
    } else {
      msg << "Found dimension size less than one in unit vector declaration"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
    }
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using stan::model::cons_list;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::nil_index_list;
using namespace stan::math;
using stan::math::pow; 
stan::math::profile_map profiles__;
static int current_statement__= 0;
static const std::vector<string> locations_array__ = {" (found before start of program)",
                                                      " (in 'string', line 9, column 2 to column 22)",
                                                      " (in 'string', line 10, column 2 to column 30)",
                                                      " (in 'string', line 13, column 2 to column 36)",
                                                      " (in 'string', line 14, column 2 to column 23)",
                                                      " (in 'string', line 16, column 5 to column 55)",
                                                      " (in 'string', line 15, column 15 to line 17, column 3)",
                                                      " (in 'string', line 15, column 2 to line 17, column 3)",
                                                      " (in 'string', line 18, column 2 to column 18)",
                                                      " (in 'string', line 21, column 2 to column 41)",
                                                      " (in 'string', line 22, column 2 to column 39)",
                                                      " (in 'string', line 2, column 2 to column 25)",
                                                      " (in 'string', line 3, column 2 to column 19)",
                                                      " (in 'string', line 4, column 2 to column 19)",
                                                      " (in 'string', line 5, column 2 to column 23)",
                                                      " (in 'string', line 6, column 20 to column 21)",
                                                      " (in 'string', line 6, column 2 to column 26)",
                                                      " (in 'string', line 10, column 20 to column 21)",
                                                      " (in 'string', line 13, column 31 to column 32)"};
#include <stan_meta_header.hpp>
class model_prior_meff final : public model_base_crtp<model_prior_meff> {
private:
  double tausq0;
  int D;
  int n;
  int sigma;
  Eigen::Matrix<double, -1, 1> s2;
 
public:
  ~model_prior_meff() { }
  
  inline std::string model_name() const final { return "model_prior_meff"; }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = %%NAME%%3 %%VERSION%%", "stancflags = "};
  }
  
  
  model_prior_meff(stan::io::var_context& context__,
                   unsigned int random_seed__ = 0,
                   std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "model_prior_meff_namespace::model_prior_meff";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 11;
      context__.validate_dims("data initialization","tausq0","double",
          context__.to_vec());
      tausq0 = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 11;
      tausq0 = context__.vals_r("tausq0")[(1 - 1)];
      current_statement__ = 11;
      current_statement__ = 11;
      check_greater_or_equal(function__, "tausq0", tausq0, 0);
      current_statement__ = 12;
      context__.validate_dims("data initialization","D","int",
          context__.to_vec());
      D = std::numeric_limits<int>::min();
      
      current_statement__ = 12;
      D = context__.vals_i("D")[(1 - 1)];
      current_statement__ = 12;
      current_statement__ = 12;
      check_greater_or_equal(function__, "D", D, 0);
      current_statement__ = 13;
      context__.validate_dims("data initialization","n","int",
          context__.to_vec());
      n = std::numeric_limits<int>::min();
      
      current_statement__ = 13;
      n = context__.vals_i("n")[(1 - 1)];
      current_statement__ = 13;
      current_statement__ = 13;
      check_greater_or_equal(function__, "n", n, 0);
      current_statement__ = 14;
      context__.validate_dims("data initialization","sigma","int",
          context__.to_vec());
      sigma = std::numeric_limits<int>::min();
      
      current_statement__ = 14;
      sigma = context__.vals_i("sigma")[(1 - 1)];
      current_statement__ = 14;
      current_statement__ = 14;
      check_greater_or_equal(function__, "sigma", sigma, 0);
      current_statement__ = 15;
      validate_non_negative_index("s2", "D", D);
      current_statement__ = 16;
      context__.validate_dims("data initialization","s2","double",
          context__.to_vec(D));
      s2 = Eigen::Matrix<double, -1, 1>(D);
      stan::math::fill(s2, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> s2_flat__;
        current_statement__ = 16;
        assign(s2_flat__, nil_index_list(), context__.vals_r("s2"),
          "assigning variable s2_flat__");
        current_statement__ = 16;
        pos__ = 1;
        current_statement__ = 16;
        for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
          current_statement__ = 16;
          assign(s2, cons_list(index_uni(sym1__), nil_index_list()),
            s2_flat__[(pos__ - 1)], "assigning variable s2");
          current_statement__ = 16;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 16;
      for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
        current_statement__ = 16;
        current_statement__ = 16;
        check_greater_or_equal(function__, "s2[sym1__]", s2[(sym1__ - 1)], 0);
      }
      current_statement__ = 17;
      validate_non_negative_index("lambda", "D", D);
      current_statement__ = 18;
      validate_non_negative_index("k", "D", D);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      num_params_r__ += 1;
      num_params_r__ += D;
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI, stan::require_vector_like_t<VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    static const char* function__ = "model_prior_meff_namespace::log_prob";
(void) function__;  // suppress unused var warning
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      local_scalar_t__ tau;
      tau = DUMMY_VAR__;
      
      current_statement__ = 1;
      tau = in__.scalar();
      current_statement__ = 1;
      if (jacobian__) {
        current_statement__ = 1;
        tau = stan::math::lb_constrain(tau, 0, lp__);
      } else {
        current_statement__ = 1;
        tau = stan::math::lb_constrain(tau, 0);
      }
      Eigen::Matrix<local_scalar_t__, -1, 1> lambda;
      lambda = Eigen::Matrix<local_scalar_t__, -1, 1>(D);
      stan::math::fill(lambda, DUMMY_VAR__);
      
      current_statement__ = 2;
      lambda = in__.vector(D);
      current_statement__ = 2;
      for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
        current_statement__ = 2;
        if (jacobian__) {
          current_statement__ = 2;
          assign(lambda, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(lambda[(sym1__ - 1)], 0, lp__),
            "assigning variable lambda");
        } else {
          current_statement__ = 2;
          assign(lambda, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(lambda[(sym1__ - 1)], 0),
            "assigning variable lambda");
        }}
      Eigen::Matrix<local_scalar_t__, -1, 1> k;
      k = Eigen::Matrix<local_scalar_t__, -1, 1>(D);
      stan::math::fill(k, DUMMY_VAR__);
      
      local_scalar_t__ meff;
      meff = DUMMY_VAR__;
      
      current_statement__ = 7;
      for (int d = 1; d <= D; ++d) {
        current_statement__ = 5;
        assign(k, cons_list(index_uni(d), nil_index_list()),
          (1 /
            (1 +
              ((((n * pow(sigma, -2)) * pow(tau, 2)) * s2[(d - 1)]) *
                pow(lambda[(d - 1)], 2)))), "assigning variable k");}
      current_statement__ = 8;
      meff = sum(subtract(1, k));
      current_statement__ = 3;
      for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
        current_statement__ = 3;
        current_statement__ = 3;
        check_greater_or_equal(function__, "k[sym1__]", k[(sym1__ - 1)], 0);}
      current_statement__ = 3;
      for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
        current_statement__ = 3;
        current_statement__ = 3;
        check_less_or_equal(function__, "k[sym1__]", k[(sym1__ - 1)], 1);}
      current_statement__ = 4;
      current_statement__ = 4;
      check_greater_or_equal(function__, "meff", meff, 0);
      {
        current_statement__ = 9;
        lp_accum__.add(cauchy_lpdf<false>(tau, 0, tausq0));
        current_statement__ = 10;
        lp_accum__.add(cauchy_lpdf<false>(lambda, 0, 1));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr>
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.resize(0);
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    static const char* function__ = "model_prior_meff_namespace::write_array";
(void) function__;  // suppress unused var warning
    (void) function__;  // suppress unused var warning
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      double tau;
      tau = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      tau = in__.scalar();
      current_statement__ = 1;
      tau = stan::math::lb_constrain(tau, 0);
      Eigen::Matrix<double, -1, 1> lambda;
      lambda = Eigen::Matrix<double, -1, 1>(D);
      stan::math::fill(lambda, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 2;
      lambda = in__.vector(D);
      current_statement__ = 2;
      for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
        current_statement__ = 2;
        assign(lambda, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_constrain(lambda[(sym1__ - 1)], 0),
          "assigning variable lambda");}
      Eigen::Matrix<double, -1, 1> k;
      k = Eigen::Matrix<double, -1, 1>(D);
      stan::math::fill(k, std::numeric_limits<double>::quiet_NaN());
      
      double meff;
      meff = std::numeric_limits<double>::quiet_NaN();
      
      vars__.emplace_back(tau);
      for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
        vars__.emplace_back(lambda[(sym1__ - 1)]);}
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      current_statement__ = 7;
      for (int d = 1; d <= D; ++d) {
        current_statement__ = 5;
        assign(k, cons_list(index_uni(d), nil_index_list()),
          (1 /
            (1 +
              ((((n * pow(sigma, -2)) * pow(tau, 2)) * s2[(d - 1)]) *
                pow(lambda[(d - 1)], 2)))), "assigning variable k");}
      current_statement__ = 8;
      meff = sum(subtract(1, k));
      current_statement__ = 3;
      for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
        current_statement__ = 3;
        current_statement__ = 3;
        check_greater_or_equal(function__, "k[sym1__]", k[(sym1__ - 1)], 0);}
      current_statement__ = 3;
      for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
        current_statement__ = 3;
        current_statement__ = 3;
        check_less_or_equal(function__, "k[sym1__]", k[(sym1__ - 1)], 1);}
      current_statement__ = 4;
      current_statement__ = 4;
      check_greater_or_equal(function__, "meff", meff, 0);
      if (emit_transformed_parameters__) {
        for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
          vars__.emplace_back(k[(sym1__ - 1)]);}
        vars__.emplace_back(meff);
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, stan::require_std_vector_t<VecVar>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void transform_inits_impl(const stan::io::var_context& context__,
                                   VecI& params_i__, VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.clear();
    vars__.reserve(num_params_r__);
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      double tau;
      tau = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      tau = context__.vals_r("tau")[(1 - 1)];
      double tau_free__;
      tau_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      tau_free__ = stan::math::lb_free(tau, 0);
      Eigen::Matrix<double, -1, 1> lambda;
      lambda = Eigen::Matrix<double, -1, 1>(D);
      stan::math::fill(lambda, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> lambda_flat__;
        current_statement__ = 2;
        assign(lambda_flat__, nil_index_list(), context__.vals_r("lambda"),
          "assigning variable lambda_flat__");
        current_statement__ = 2;
        pos__ = 1;
        current_statement__ = 2;
        for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
          current_statement__ = 2;
          assign(lambda, cons_list(index_uni(sym1__), nil_index_list()),
            lambda_flat__[(pos__ - 1)], "assigning variable lambda");
          current_statement__ = 2;
          pos__ = (pos__ + 1);}
      }
      Eigen::Matrix<double, -1, 1> lambda_free__;
      lambda_free__ = Eigen::Matrix<double, -1, 1>(D);
      stan::math::fill(lambda_free__, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 2;
      for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
        current_statement__ = 2;
        assign(lambda_free__, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_free(lambda[(sym1__ - 1)], 0),
          "assigning variable lambda_free__");}
      vars__.emplace_back(tau_free__);
      for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
        vars__.emplace_back(lambda_free__[(sym1__ - 1)]);}
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__.clear();
    names__.emplace_back("tau");
    names__.emplace_back("lambda");
    names__.emplace_back("k");
    names__.emplace_back("meff");
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.clear();
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(D)});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(D)});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "tau");
    for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "lambda" + '.' + std::to_string(sym1__));
      }}
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "k" + '.' + std::to_string(sym1__));
        }}
      param_names__.emplace_back(std::string() + "meff");
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "tau");
    for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "lambda" + '.' + std::to_string(sym1__));
      }}
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= D; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "k" + '.' + std::to_string(sym1__));
        }}
      param_names__.emplace_back(std::string() + "meff");
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"tau\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"lambda\",\"type\":{\"name\":\"vector\",\"length\":" << D << "},\"block\":\"parameters\"},{\"name\":\"k\",\"type\":{\"name\":\"vector\",\"length\":" << D << "},\"block\":\"transformed_parameters\"},{\"name\":\"meff\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"tau\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"lambda\",\"type\":{\"name\":\"vector\",\"length\":" << D << "},\"block\":\"parameters\"},{\"name\":\"k\",\"type\":{\"name\":\"vector\",\"length\":" << D << "},\"block\":\"transformed_parameters\"},{\"name\":\"meff\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"}]";
    return s__.str();
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      std::vector<double> vars_vec(vars.size());
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i) {
        vars.coeffRef(i) = vars_vec[i];
      }
    }
    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      write_array_impl(base_rng, params_r, params_i, vars, emit_transformed_parameters, emit_generated_quantities, pstream);
    }
    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
  
    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits_impl(context, params_i, params_r_vec, pstream);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i) {
        params_r.coeffRef(i) = params_r_vec[i];
      }
    }
    inline void transform_inits(const stan::io::var_context& context,
                                std::vector<int>& params_i,
                                std::vector<double>& vars,
                                std::ostream* pstream = nullptr) const final {
      transform_inits_impl(context, params_i, vars, pstream);
    }        
};
}
using stan_model = model_prior_meff_namespace::model_prior_meff;
#ifndef USING_R
// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_prior_meff_namespace::profiles__;
}
#endif
#endif
