#include "smcpp.hpp"

#include "mpreal.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/MatrixFunctions>

#include <cereal/archives/portable_binary.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/map.hpp>

#include <gmpxx.h>

#include <fcntl.h> // for open()

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#ifdef NO_CHECK_NAN
#define CHECK_NAN(x)
#define CHECK_NAN_OR_NEGATIVE(x)
#else
#define CHECK_NAN(x) check_nan(x, __FILE__, __LINE__)
#define CHECK_NAN_OR_NEGATIVE(x) check_nan(x, __FILE__, __LINE__); check_negative(x, __FILE__, __LINE__);
#endif

#define ERROR std::cout << "smcpp: "
#define DEBUG1 std::cout << "smcpp: "


using adouble_base_type = double;
using adouble_t = Eigen::Matrix<adouble_base_type, Eigen::Dynamic, 1>;
using adouble = Eigen::AutoDiffScalar<adouble_t>;
using ParameterVector = std::vector<std::vector<adouble>>;


void check_nan(const double x, const char* file, const int line);
void check_nan(const adouble &x, const char* file, const int line);

void check_negative(const adouble x, const char* file, const int line);
void check_negative(const double x, const char* file, const int line);


template <typename T>
inline T doubly_compensated_summation(const std::vector<T> &x)
{
  if (x.size() == 0)
    return 0.0;
  T s = x[0];
  T c = 0.0;
  T y, u, v, t, z;
  for (unsigned int i = 1; i < x.size(); ++i)
  {
    y = c + x[i];
    u = x[i] - (y - c);
    t = y + s;
    v = y - (t - s);
    z = u + v;
    s = t + z;
    c = z - (s - t);
  }
  return s;
}


inline double toDouble(const adouble &a) { return a.value(); }
inline double toDouble(const double &d) { return d; }


namespace Eigen {
// Allow for casting of adouble matrices to double
namespace internal
{
template <>
struct cast_impl<adouble, float>
{
  static inline double run(const adouble &x)
  {
    return x.value();
  }
};
template <>
struct cast_impl<adouble, double>
{
  static inline double run(const adouble &x)
  {
    return x.value();
  }
};
}

// Copied from Eigen's AutoDiffScalar.h
#define EIGEN_AUTODIFF_DECLARE_GLOBAL_UNARY(FUNC,CODE) \
  template<typename DerType> \
  inline const Eigen::AutoDiffScalar< \
  EIGEN_EXPR_BINARYOP_SCALAR_RETURN_TYPE(typename Eigen::internal::remove_all<DerType>::type, typename Eigen::internal::traits<typename Eigen::internal::remove_all<DerType>::type>::Scalar, product) > \
  FUNC(const Eigen::AutoDiffScalar<DerType>& x) { \
    using namespace Eigen; \
    EIGEN_UNUSED typedef typename Eigen::internal::traits<typename Eigen::internal::remove_all<DerType>::type>::Scalar Scalar; \
    CODE; \
  }

using std::log1p;
using std::cosh;
using std::sinh;
using mpfr::exp;
using mpfr::sinh;
using mpfr::cosh;

EIGEN_AUTODIFF_DECLARE_GLOBAL_UNARY(expm1,
                                    using std::exp;
                                        using std::expm1;
                                        Scalar expm1x = expm1(x.value());
                                        Scalar expx = exp(x.value());
                                        return Eigen::MakeAutoDiffScalar(expm1x, x.derivatives() * expx);
)

EIGEN_AUTODIFF_DECLARE_GLOBAL_UNARY(log1p,
                                    using std::log1p;
                                        Scalar log1px = log1p(x.value());
                                        return Eigen::MakeAutoDiffScalar(log1px, x.derivatives() * (Scalar(1) / (Scalar(1) + x.value())));
)

#undef EIGEN_AUTODIFF_DECLARE_GLOBAL_UNARY

};


template <typename Derived>
void check_nan(const Eigen::DenseBase<Derived> &M, const char* file, const int line)
{
  for (int i = 0; i < M.rows(); ++i)
    for (int j = 0; j < M.cols(); ++j)
    {
      try
      {
        check_nan(M.coeff(i, j), file, line);
      }
      catch (std::runtime_error)
      {
//        CRITICAL << "rows:" << M.rows() << " cols:" << M.cols()
//                 << " M(" << i << "," << j << ")=" << M.coeff(i, j);
        throw;
      }
    }
}


template <typename T>
class PiecewiseConstantRateFunction
{
public:
  PiecewiseConstantRateFunction(const std::vector<std::vector<adouble>>, const std::vector<double>);
  PiecewiseConstantRateFunction(const PiecewiseConstantRateFunction &other) :
      PiecewiseConstantRateFunction(other.params, other.hidden_states) {}
  T zero() const;
  T R(const T) const;
  T Rinv(const T) const;
  T R_integral(const double, const double, const T) const;
  std::vector<T> average_coal_times() const;
  T random_time(const double, const double, const long long) const;
  T random_time(const double, const double, const double, std::mt19937 &) const;
  int getNder() const { return nder; }
  void print_debug() const;

  void tjj_double_integral_above(const int, long, std::vector<Matrix<T> > &) const;
  void tjj_double_integral_below(const int, const int, Matrix<T>&) const;

  // Getters
  const std::vector<double>& getHiddenStates() const { return hidden_states; }
  const std::vector<double>& getTs() const { return ts; }
  const std::vector<int>& getHsIndices() const { return hs_indices; }
  const std::vector<T>& getRrng() const { return Rrng; }
  const std::vector<T>& getAda() const { return ada; }

  friend std::ostream& operator<<(std::ostream& os, const PiecewiseConstantRateFunction& pexp)
  {
//    os << pexp.ts << std::endl;
//    os << pexp.ada << std::endl;
    return os;
  }

private:
  const std::vector<std::vector<adouble>> params;
  const int nder;
  int K;
  std::vector<T> ada;
  std::vector<double> s;
  std::vector<double> ts;
  std::vector<T> Rrng;
  const std::vector<double> hidden_states;
  std::vector<int> hs_indices;
  // Methods
  void compute_antiderivative();
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// mpq_support.h ///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace Eigen {
template<>
struct NumTraits<mpq_class> : NumTraits<long long> // permits to get the epsilon, dummy_precision, lowest, highest functions
{
  typedef mpq_class Real;
  typedef mpq_class Nested;
  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 1,
    ReadCost = 1,
    AddCost = 3,
    MulCost = 3
  };
};
namespace internal
{
template <>
struct scalar_product_traits<double, mpq_class>
{
  typedef double ReturnType;
};
template <>
struct scalar_product_traits<adouble, mpq_class>
{
  typedef adouble ReturnType;
};
template <>
struct cast_impl<mpq_class, adouble>
{
  static inline adouble run(const mpq_class &x)
  {
    return adouble(mpq_get_d(x.get_mpq_t()));
  }
};
template <>
struct cast_impl<mpq_class, double>
{
  static inline double run(const mpq_class &x)
  {
    return mpq_get_d(x.get_mpq_t());
  }
};
}
}

typedef Eigen::Matrix<mpq_class, Eigen::Dynamic, Eigen::Dynamic> MatrixXq;
typedef Eigen::Matrix<mpq_class, Eigen::Dynamic, 1> VectorXq;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// moran_eigensystem.h /////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct MoranEigensystem
{
  MoranEigensystem(const int n) : U(n + 1, n + 1), Uinv(n + 1, n + 1), D(n + 1)
  {
    U.setZero();
    Uinv.setZero();
    D.setZero();
  }
  MatrixXq U, Uinv;
  VectorXq D;
};

Eigen::SparseMatrix<mpq_class, Eigen::RowMajor> moran_rate_matrix(int);
Eigen::SparseMatrix<mpq_class, Eigen::RowMajor> modified_moran_rate_matrix(int, int, int);
MoranEigensystem& compute_moran_eigensystem(int);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// matrix_cache.h //////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct MatrixCache
{
  template <class Archive>
  void serialize(Archive & ar)
  {
    ar(X0, X2, M0, M1);
  }
  Matrix<double> X0, X2, M0, M1;
};

MatrixCache& cached_matrices(int);



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








template <typename T>
class ConditionedSFS
{
public:
  virtual ~ConditionedSFS() = default;
  virtual std::vector<Matrix<T> > compute(const PiecewiseConstantRateFunction<T> &) = 0;
};

template <typename T>
class OnePopConditionedSFS : public ConditionedSFS<T>
{
public:
  OnePopConditionedSFS(int);
  std::vector<Matrix<T> > compute(const PiecewiseConstantRateFunction<T> &);

  std::vector<Matrix<T> > compute_below(const PiecewiseConstantRateFunction<T> &) const;
  std::vector<Matrix<T> > compute_above(const PiecewiseConstantRateFunction<T> &) const;

private:
  const int n;
  const MoranEigensystem mei;
  const MatrixCache mcache;
  const Matrix<double> Uinv_mp0, Uinv_mp2;
};



















void check_nan(const double x, const char* file, const int line)
{
  std::string s;
  if (std::isnan(x))
    s = "nan";
  else if (std::isinf(x))
    s = "inf";
  else
    return;
  s += " detected at ";
  s += file;
  s += ":";
  s += std::to_string(line);
//#pragma omp critical(stacktrace)
//  {
//    CRITICAL << s;
//    print_stacktrace();
//  }
  throw std::runtime_error(s);
}

void check_nan(const adouble &x, const char* file, const int line)
{
  check_nan(x.value(), file, line);
  check_nan(x.derivatives(), file, line);
}


void check_negative(const adouble x, const char* file, const int line) {
  check_negative(x.value(), file, line);
}

void check_negative(const double x, const char* file, const int line)
{
  if (x > -1e-16)
    return;
  std::string s = "negative value detected at ";
  s += file;
  s += ":";
  s += std::to_string(line);
//#pragma omp critical(stacktrace)
//  {
//    CRITICAL << s;
//    print_stacktrace();
//  }
  throw std::runtime_error(s);
}


constexpr long nC2(int n) { return n * (n - 1) / 2; }

template <typename T>
inline T _conv(const adouble x);

template <typename T>
inline std::vector<T> _vconv(const std::vector<adouble> v);

template <>
inline double _conv(const adouble x) { return x.value(); }

template <>
inline adouble _conv(const adouble x) { return x; }

template <>
inline std::vector<adouble> _vconv(const std::vector<adouble> v) { return v; }

template <>
inline std::vector<double> _vconv(const std::vector<adouble> v)
{
  std::vector<double> ret;
  for (adouble x : v)
    ret.push_back(x.value());
  return ret;
}

template <typename T>
PiecewiseConstantRateFunction<T>::PiecewiseConstantRateFunction(
    const std::vector<std::vector<adouble>> params,
    const std::vector<double> hidden_states) :
    params(params),
    nder(params[0][0].derivatives().size()),
    K(params[0].size()),
    ada(_vconv<T>(params[0])),
    s(_vconv<double>(params[1])),
    ts(K + 1), Rrng(K),
    hidden_states(hidden_states)
{
  for (auto &pp : params)
    if (pp.size() != params[0].size())
      throw std::runtime_error("all params must have same size");
  // Final piece is required to be flat.
  ts[0] = 0.;
  Rrng[0] = 0.;
  // These constant values need to have compatible derivative shape
  // with the calculated values.
  // Fix last piece to be constant
  for (int k = 0; k < K; ++k)
  {
    ada[k] = 1. / ada[k];
    ts[k + 1] = ts[k] + s[k];
  }
  ts[K] = INFINITY;

  for (double h : hidden_states)
  {
    if (std::isinf(h))
    {
      hs_indices.push_back(ts.size() - 1);
      continue;
    }
    std::vector<double>::iterator ti = std::upper_bound(ts.begin(), ts.end(), h) - 1;
    int ip = ti - ts.begin();
    if (std::abs(*ti - h) < 1e-8)
      // if (ts[ip] == h)
      hs_indices.push_back(ip);
    else if (ti + 1 < ts.end() and std::abs(*(ti + 1) - h) < 1e-8)
      hs_indices.push_back(ip + 1);
    else
    {
      ts.insert(ti + 1, h);
      ada.insert(ada.begin() + ip + 1, ada[ip]);
      CHECK_NAN(ada[ip + 1]);
      CHECK_NAN(ts[ip + 1]);
      hs_indices.push_back(ip + 1);
    }
  }
  K = ada.size();
  Rrng.resize(K + 1);
  compute_antiderivative();
}

template <typename T>
inline T _double_integral_below_helper(const int rate, const double tsm, const double tsm1, const T ada,
                                       const T Rrng, const T log_denom)
{
  if (ada == 0)
    return 0.;
  const int l1r = 1 + rate;
  const double l1rinv = 1. / (double)l1r;
  double diff = tsm1 - tsm;
  T adadiff = ada * diff;
  if (rate == 0)
  {
    if (tsm1 == INFINITY)
      return exp(-Rrng - log_denom) / ada;
    else
      return exp(-Rrng - log_denom) * (1. - exp(-adadiff) * (1. + adadiff)) / ada;
  }
  if (tsm1 == INFINITY)
    return exp(-l1r * Rrng - log_denom) * (1. - l1rinv) / (rate * ada);
  return exp(-l1r * Rrng - log_denom) * (expm1(-l1r * adadiff) * l1rinv - expm1(-adadiff)) / (rate * ada);
}

template <typename T>
inline T _double_integral_above_helper(const int rate, const int lam, const double tsm,
                                       const double tsm1, const T ada, const T Rrng, const T log_coef)
{
  if (ada == 0)
    return 0.;
  double diff = tsm1 - tsm;
  T adadiff = ada * diff;
  long l1 = lam + 1;
  if (rate == 0)
  {
    return exp(-l1 * Rrng + log_coef) * (expm1(-l1 * adadiff) + l1 * adadiff) / l1 / l1 / ada;
  }
  if (l1 == rate)
  {
    if (tsm1 == INFINITY)
      return exp(-rate * Rrng + log_coef) / rate / rate / ada;
    return exp(-rate * Rrng + log_coef) * (1 - exp(-rate * adadiff) * (1 + rate * adadiff)) / rate / rate / ada;
  }
  if (tsm1 == INFINITY)
    return exp(-l1 * Rrng + log_coef) / l1 / rate / ada;
  // return -exp(-l1 * _Rrng + log_coef) * (expm1(-l1 * adadiff) / l1 + (exp(-rate * adadiff) - exp(-l1 * adadiff)) / (l1 - rate)) / rate / _ada;
  if (rate < l1)
    return -exp(-l1 * Rrng + log_coef) * (expm1(-l1 * adadiff) / l1 + (exp(-rate * adadiff) * -expm1(-(l1 - rate) * adadiff) / (l1 - rate))) / rate / ada;
  else
    return -exp(-l1 * Rrng + log_coef) *
           (
               expm1(-l1 * adadiff) / l1 +
               (exp(-l1 * adadiff) * expm1(-(rate - l1) * adadiff) / (l1 - rate))
           ) / rate / ada;
}

template <typename T>
void PiecewiseConstantRateFunction<T>::print_debug() const
{
  std::vector<std::pair<std::string, std::vector<T> > > arys = {{"ada", ada}, {"Rrng", Rrng}};
  for (auto p : arys)
  {
//    CRITICAL << p.first << "\n";
    for (adouble x : p.second)
    {
//      CRITICAL << x.value()
//               << "::" << x.derivatives().transpose()
//               << "\n";
    }
  }
}

template <typename T>
void PiecewiseConstantRateFunction<T>::compute_antiderivative()
{
  Rrng[0] = 0.;
  for (int k = 0; k < K; ++k)
    Rrng[k + 1] = Rrng[k] + ada[k] * (ts[k + 1] - ts[k]);
}

template <typename T>
T PiecewiseConstantRateFunction<T>::R_integral(const double a,
                                               const double b, const T log_denom) const
{
  // int_a^b exp(-R(t)) dt
  int ip_a = std::upper_bound(ts.begin(), ts.end(), a) - 1 - ts.begin();
  int ip_b = std::upper_bound(ts.begin(), ts.end(), b) - 1 - ts.begin();
  // If b == inf the comparison is always false so a special case is needed.
  ip_b = std::isinf(b) ? ts.size() - 2 : ip_b;
  T ret = 0., r, Rleft;
  double left, right, diff;
  for (int i = ip_a; i < ip_b + 1; ++i)
  {
    left = std::max(a, ts[i]);
    right = std::min(b, ts[i + 1]);
    diff = right - left;
    Rleft = R(left);
    r = exp(-(Rleft + log_denom));
    if (ada[i] > 0.)
    {
      if (!std::isinf(toDouble(diff)))
        r *= -expm1(-diff * ada[i]);
      r /= ada[i];
    }
    else
      r *= diff;
    CHECK_NAN_OR_NEGATIVE(r);
    ret += r;
  }
  return ret;
}


template <typename T>
inline T _single_integral(const int rate, const double tsm, const double tsm1,
                          const T ada, const T Rrng, const T log_coef)
{
  // = int_ts[m]^ts[m+1] exp(-rate * R(t)) dt
  const int c = rate;
  if (rate == 0)
    return exp(log_coef) * (tsm1 - tsm);
  T ret = exp(-c * Rrng + log_coef);
  if (tsm1 < INFINITY)
    ret *= -expm1(-c * ada * (tsm1 - tsm));
  ret /= ada * c;
  CHECK_NAN_OR_NEGATIVE(ret);
  return ret;
}

template <typename T>
void PiecewiseConstantRateFunction<T>::tjj_double_integral_above(
    const int n, long jj, std::vector<Matrix<T> > &C) const
{
  T tmp;
  const T z = zero();
  long lam = nC2(jj) - 1;
  // Now calculate with hidden state integration limits
  for (unsigned int h = 0; h < hs_indices.size() - 1; ++h)
  {
    C[h].row(jj - 2).fill(z);
    T Rh = Rrng[hs_indices[h]];
    T Rh1 = Rrng[hs_indices[h + 1]];
    T log_denom = -Rh;
    if (Rh1 != INFINITY)
      log_denom += log(-expm1(-(Rh1 - Rh)));
    for (int m = hs_indices[h]; m < hs_indices[h + 1]; ++m)
    {
      for (int j = 2; j < n + 2; ++j)
      {
        long rate = nC2(j);
        tmp = _double_integral_above_helper<T>(rate, lam, ts[m], ts[m + 1], ada[m], Rrng[m], -log_denom);
        try
        {
          CHECK_NAN_OR_NEGATIVE(tmp);
          CHECK_NAN_OR_NEGATIVE(C[h](jj - 2, j - 2));
        }
        catch (std::runtime_error)
        {
//          CRITICAL << "nan detected:\n j=" << j << "m=" << m << " rate=" << rate
//                   << " lam=" << lam << " ts[m]=" << ts[m] << " ts[m + 1]="
//                   << ts[m + 1] << " ada[m]=" << ada[m] << " Rrng[m]=" << Rrng[m]
//                   << " log_denom=" << log_denom;
//          CRITICAL << "tmp=" << tmp;
//          CRITICAL << "C[h](jj - 2, j - 2)= " << C[h](jj - 2, j - 2);
//          CRITICAL << "h=" << h;
//          CRITICAL << "I am eta: ";
          print_debug();
          throw;
        }
        C[h](jj - 2, j - 2) += tmp;
        T log_coef = -log_denom, fac;
        long rp = lam + 1 - rate;
        T Rm1 = Rrng[m + 1];
        T Rm = Rrng[m];
        if (rp == 0)
          fac = Rm1 - Rm;
        else
        {
          if (rp < 0)
          {
            if (-rp * (Rm1 - Rm) > 20)
            {
              log_coef += -rp * Rm1;
              fac = -1. / rp;
            }
            else
            {
              log_coef += -rp * Rm;
              fac = -expm1(-rp * (Rm1 - Rm)) / rp;
            }
          }
          else
          {
            if (-rp * (Rm - Rm1) > 20)
            {
              log_coef += -rp * Rm;
              fac = 1. / rp;
            }
            else
            {
              log_coef += -rp * Rm1;
              fac = expm1(-rp * (Rm - Rm1)) / rp;
            }
          }
        }
        for (int k = m + 1; k < K; ++k)
        {
          T si = _single_integral(rate, ts[k], ts[k + 1], ada[k], Rrng[k], log_coef) * fac;
          C[h](jj - 2, j - 2) += si;
          CHECK_NAN_OR_NEGATIVE(C[h](jj - 2, j - 2));
        }
        CHECK_NAN_OR_NEGATIVE(C[h](jj - 2, j - 2));
      }
    }
  }
}

template <typename T>
void PiecewiseConstantRateFunction<T>::tjj_double_integral_below(
    const int n, const int h, Matrix<T> &tgt) const
{
//  DEBUG1 << "in tjj_double_integral_below";
  T Rh = Rrng[hs_indices[h]];
  T Rh1 = Rrng[hs_indices[h + 1]];
  T log_denom = -Rh;
  if (Rh1 != INFINITY)
    log_denom += log(-expm1(-(Rh1 - Rh)));
  for (int m = hs_indices[h]; m < hs_indices[h + 1]; ++m)
  {
    T Rm = Rrng[m];
    T Rm1 = Rrng[m + 1];
    Vector<T> ts_integrals(n + 1);
    T log_coef = -Rm;
    T fac = 1.;
    if (m < K - 1)
      fac = -expm1(-(Rm1 - Rm));
    for (int j = 2; j < n + 3; ++j)
    {
      long rate = nC2(j) - 1;
      ts_integrals(j - 2) = _double_integral_below_helper<T>(rate, ts[m], ts[m + 1], ada[m], Rrng[m], log_denom);
      for (int k = 0; k < m; ++k)
      {
        T _c = log_coef - log_denom;
        ts_integrals(j - 2) += fac * _single_integral(rate, ts[k], ts[k + 1], ada[k], Rrng[k], _c);
      }
      CHECK_NAN_OR_NEGATIVE(ts_integrals(j - 2));
    }
    tgt.row(h) += ts_integrals.transpose();
  }
//  DEBUG1 << "exiting tjj_double_integral_below";
}

template <typename T>
T exp1_conditional(T a, T b, std::mt19937 &gen)
{
  // If X ~ Exp(1),
  // P(X < x | a <= X <= b) = (e^-a - e^-x) / (e^-a - e^-b)
  // so P^-1(y) = -log(e^-a - (e^-a - e^-b) * y)
  //            = -log(e^-a(1 - (1 - e^-(b-a)) * y)
  //            = a - log(1 - (1 - e^-(b-a)) * y)
  //            = a - log(1 + expm1(-(b-a)) * y)
  double unif = std::uniform_real_distribution<double>{0.0, 1.0}(gen);
  if (std::isinf(toDouble(b)))
    return a - log1p(-unif);
  else
    return a - log1p(expm1(-(b - a)) * unif);
}

template <typename T>
T PiecewiseConstantRateFunction<T>::random_time(const double a, const double b, const long long seed) const
{
  std::mt19937 gen(seed);
  return random_time(1., a, b, gen);
}

template <typename T>
T PiecewiseConstantRateFunction<T>::random_time(const double fac, const double a, const double b, std::mt19937 &gen) const
{
  T Rb;
  if (b == INFINITY)
    Rb = INFINITY;
  else
    Rb = R(b);
  return Rinv(exp1_conditional(R(a), Rb, gen) / fac);
}


template <typename T>
std::vector<T> PiecewiseConstantRateFunction<T>::average_coal_times() const
{
  std::vector<T> ret;
  for (unsigned int i = 1; i < hidden_states.size(); ++i)
  {
    // discretize by expected coalescent time within each hidden state
    // e_coal = \int_t0^t1 t eta(t) exp(-R(t))
    //        = t0 e^(-R(t0)) - t1 e^(-R(t1)) + \int
    if (Rrng[hs_indices[i - 1]] == Rrng[hs_indices[i]])
    {
      // Handle the case of infinite population size => no chance
      // of coalescence. This can happen in a two population model
      // where the distinguished lineages come from separate
      // populations.
      ret.push_back(std::numeric_limits<double>::quiet_NaN());
      continue;
    }
    T log_denom = -Rrng[hs_indices[i - 1]];
    bool inf = std::isinf(toDouble(ts[hs_indices[i]]));
    if (!inf)
      log_denom += log(-expm1(-(Rrng[hs_indices[i]] - Rrng[hs_indices[i - 1]])));
    T x = hidden_states[i - 1] * exp(-((Rrng[hs_indices[i - 1]]) + log_denom)) +
          R_integral(ts[hs_indices[i - 1]], ts[hs_indices[i]], log_denom);
    if (!inf)
      x -= hidden_states[i] * exp(-((Rrng[hs_indices[i]]) + log_denom));
    CHECK_NAN(x);
    ret.push_back(x);
    if (ret.back() > hidden_states[i] or ret.back() < hidden_states[i - 1])
      throw std::runtime_error("erroneous average coalescence time");
  }
  return ret;
}

template <typename T>
T PiecewiseConstantRateFunction<T>::R(const T t) const
{
  std::vector<double>::const_iterator ti = std::upper_bound(ts.begin(), ts.end(), toDouble(t)) - 1;
  int ip = ti - ts.begin();
  return Rrng[ip] + ada[ip] * (t - *ti);
}


template <typename T>
T PiecewiseConstantRateFunction<T>::Rinv(const T y) const
{
  typename std::vector<T>::const_iterator R = std::upper_bound(Rrng.begin(), Rrng.end(), y) - 1;
  int ip = R - Rrng.begin();
  return (y - *R) / ada[ip] + ts[ip];
}

template <>
double PiecewiseConstantRateFunction<double>::zero() const
{
  return 0.;
}

template <>
adouble PiecewiseConstantRateFunction<adouble>::zero() const
{
  return adouble(0., Vector<double>::Zero(nder));
}

template class PiecewiseConstantRateFunction<double>;
template class PiecewiseConstantRateFunction<adouble>;

const mpq_class mpq_1(2, 1);
const mpq_class mpq_2(2, 1);




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// moran_eigensystem.cpp////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::SparseMatrix<mpq_class, Eigen::RowMajor> moran_rate_matrix(int N)
{
  Eigen::SparseMatrix<mpq_class, Eigen::RowMajor> ret(N + 1, N + 1);
  for (int i = 0; i < N + 1; ++i)
  {
    mpq_class sm = 0, b;
    if (i > 0)
    {
      b = i * (N - i) / mpq_2;
      ret.insert(i, i - 1) = b;
      sm += b;
    }
    if (i < N)
    {
      b = i * (N - i) / mpq_2;
      ret.insert(i, i + 1) = b;
      sm += b;
    }
    ret.insert(i, i) = -sm;
  }
  return ret;
}

Eigen::SparseMatrix<mpq_class, Eigen::RowMajor> modified_moran_rate_matrix(int N, int a, int na)
{
  Eigen::SparseMatrix<mpq_class, Eigen::RowMajor> ret(N + 1, N + 1);
  for (int i = 0; i < N + 1; ++i)
  {
    mpq_class sm = 0, b;
    if (i > 0)
    {
      b = (na - a) * i + i * (N - i) / mpq_2;
      ret.insert(i, i - 1) = b;
      sm += b;
    }
    if (i < N)
    {
      b = a * (N - i) + i * (N - i) / mpq_2;
      ret.insert(i, i + 1) = b;
      sm += b;
    }
    ret.insert(i, i) = -sm;
  }
  return ret;
}

VectorXq solve(const Eigen::SparseMatrix<mpq_class, Eigen::RowMajor> &M)
// VectorXq solve(const Eigen::SparseMatrixBase<Derived> &M)
{
  int n = M.rows();
  VectorXq ret(n);
  ret.setZero();
  ret(n - 1) = 1;
  for (int i = n - 2; i > -1; --i)
    ret(i) = (M.row(i + 1) * ret).sum() / -M.coeff(i + 1, i);
  return ret;
}

// This function is not thread safe. Do not call from multiple threads.
std::map<int, MoranEigensystem> _memo;
MoranEigensystem& compute_moran_eigensystem(int n)
{
  if (_memo.count(n) == 0)
  {
    Eigen::SparseMatrix<mpq_class, Eigen::RowMajor> M = modified_moran_rate_matrix(n, 0, 2),
        Mt, I(n + 1, n + 1), A;
    Mt = M.transpose();
    MoranEigensystem ret(n);
    ret.Uinv(0, 0) = mpq_1;
    I.setIdentity();
    for (int k = 2; k < n + 3; ++k)
    {
      mpq_class rate = -(k * (k - 1) / 2 - 1);
      ret.D(k - 2) = rate;
      A = M - rate * I;
      ret.U.col(k - 2) = solve(A);
      if (k > 2)
      {
        A = Mt - rate * I;
        ret.Uinv.row(k - 2).tail(n) = solve(A.bottomRightCorner(n, n));
        ret.Uinv(k - 2, 0) = -ret.Uinv(k - 2, 1) * A.coeff(0, 1) / A.coeff(0, 0);
      }
    }
    VectorXq D1 = (ret.Uinv * ret.U).diagonal().cwiseInverse();
    ret.U = ret.U * D1.asDiagonal();
    _memo.emplace(n, ret);
  }
  return _memo.at(n);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// matrix_cache.cpp ////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace Eigen
{
template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> inline
void
save(Archive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> const & m)
{
  int32_t rows = m.rows();
  int32_t cols = m.cols();
  ar(rows);
  ar(cols);
  ar(cereal::binary_data(m.data(), rows * cols * sizeof(_Scalar)));
}

template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> inline
void
load(Archive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & m)
{
  int32_t rows;
  int32_t cols;
  ar(rows);
  ar(cols);

  m.resize(rows, cols);

  ar(cereal::binary_data(m.data(), static_cast<std::size_t>(rows * cols * sizeof(_Scalar))));
}
}

static std::map<int, MatrixCache> cache;
static std::string store_location;

void smcpp_init_cache()
{
  store_location = (std::filesystem::temp_directory_path() / "smcpp_matrices.dat").string();
  int fd;
  const std::string lockfile = store_location + ".lock";
  fd = open(lockfile.c_str(), O_WRONLY | O_CREAT, 0600);
  struct flock oflock;
  oflock.l_type = F_WRLCK;/*Write lock*/
  oflock.l_whence = SEEK_SET;
  oflock.l_start = 0;
  oflock.l_len = 0;/*Lock whole file*/
  if (fd == -1 || fcntl(fd, F_SETLK, &oflock) == -1)
  {
//    DEBUG1 << "Couldn't acquire lock in init_cache()";
    return;
  }
//  DEBUG1 << "Acquired lock in init_cache()";
  std::ifstream in(store_location, std::ios::binary);
  if (in)
  {
    cereal::PortableBinaryInputArchive iarchive(in);
    iarchive(cache);
  }
  oflock.l_type = F_UNLCK;
  if (fcntl(fd, F_UNLCK, &oflock) == -1)
    ERROR << "Couldn't release lock in store_cache()";
  close(fd);
  remove(lockfile.c_str());
  DEBUG1 << "init_cache() successful\n";
}

void store_cache()
{
  DEBUG1 << "storing cache: " << store_location << '\n';
  int fd;
  const std::string lockfile = store_location + ".lock";
  fd = open(lockfile.c_str(), O_WRONLY | O_CREAT, 0600);
  struct flock oflock;
  oflock.l_type = F_WRLCK; /*Write lock*/
  oflock.l_whence = SEEK_SET;
  oflock.l_start = 0;
  oflock.l_len = 0;/*Lock whole file*/
  if (fd == -1 || fcntl(fd, F_SETLK, &oflock) == -1)
  {
    ERROR << "Couldn't acquire lock in store_cache()";
    return;
  }
//  DEBUG1 << "Acquired lock in store_cache()";
  std::ofstream out(store_location, std::ios::binary);
  if (out)
  {
    cereal::PortableBinaryOutputArchive oarchive(out);
    oarchive(cache);
  }
  else
  {
    ERROR << "could not open cache file for storage";
  }
  oflock.l_type = F_UNLCK;
  if (fcntl(fd, F_UNLCK, &oflock) == -1)
    ERROR << "Couldn't release lock in store_cache()";
  close(fd);
  remove(lockfile.c_str());
  DEBUG1 << "store_cache() successful\n";
}

typedef struct { MatrixXq coeffs; } below_coeff;

std::map<int, below_coeff> below_coeffs_memo;
below_coeff compute_below_coeffs(int n)
{
  if (below_coeffs_memo.count(n) == 0)
  {
//    DEBUG1 << "Computing below_coeffs";
    below_coeff ret;
    MatrixXq mlast;
    for (int nn = 2; nn < n + 3; ++nn)
    {
      MatrixXq mnew(n + 1, nn - 1);
      mnew.col(nn - 2).setZero();
      mnew(nn - 2, nn - 2) = 1;
#pragma omp parallel for
      for (int k = nn - 1; k > 1; --k)
      {
        long denom = (nn + 1) * (nn - 2) - (k + 1) * (k - 2);
        mpq_class c1((nn + 1) * (nn - 2), denom);
        mnew.col(k - 2) = mlast.col(k - 2) * c1;
      }
      for (int k = nn - 1; k > 1; --k)
      {
        long denom = (nn + 1) * (nn - 2) - (k + 1) * (k - 2);
        mpq_class c2((k + 2) * (k - 1), denom);
        mnew.col(k - 2) -= mnew.col(k - 1) * c2;
      }
      mlast = mnew;
    }
    ret.coeffs = mlast;
    below_coeffs_memo.emplace(n, ret);
  }
  return below_coeffs_memo.at(n);
}

std::map<std::array<int, 3>, mpq_class> _Wnbj_memo;
mpq_class calculate_Wnbj(int n, int b, int j)
{
  switch (j)
  {
  case 2:
    return mpq_class(6, n + 1);
  case 3:
    if (n == 2 * b) return 0;
    return mpq_class(30 * (n - 2 * b), (n + 1) * (n + 2));
  default:
    std::array<int, 3> key = {n, b, j};
    if (_Wnbj_memo.count(key) == 0)
    {
      int jj = j - 2;
      mpq_class c1(-(1 + jj) * (3 + 2 * jj) * (n - jj), jj * (2 * jj - 1) * (n + jj + 1));
      mpq_class c2((3 + 2 * jj) * (n - 2 * b), jj * (n + jj + 1));
      mpq_class ret = calculate_Wnbj(n, b, jj) * c1;
      ret += calculate_Wnbj(n, b, jj + 1) * c2;
      _Wnbj_memo[key] = ret;
    }
    return _Wnbj_memo[key];
  }
}

std::map<std::array<int, 3>, mpq_class> pnkb_dist_memo;
mpq_class pnkb_dist(int n, int m, int l1)
{
  // Probability that lineage 1 has size |L_1|=l1 below tau,
  // the time at which 1 and 2 coalesce, when there are k
  // undistinguished lineages remaining, in a sample of n
  // undistinguished (+2 distinguished) lineages overall.
  std::array<int, 3> key{n, m, l1};
  if (pnkb_dist_memo.count(key) == 0)
  {
    mpz_class binom1, binom2;
    mpz_bin_uiui(binom1.get_mpz_t(), n + 2 - l1, m + 1);
    mpz_bin_uiui(binom2.get_mpz_t(), n + 3, m + 3);
    mpq_class ret(binom1, binom2);
    ret *= l1;
    pnkb_dist_memo[key] = ret;
  }
  return pnkb_dist_memo[key];
}

std::map<std::array<int, 3>, mpq_class> pnkb_undist_memo;
mpq_class pnkb_undist(int n, int m, int l3)
{
  // Probability that undistinguished lineage has size |L_1|=l1 below tau,
  // the time at which 1 and 2 coalesce, when there are k
  // undistinguished lineages remaining, in a sample of n
  // undistinguished (+2 distinguished) lineages overall.
  std::array<int, 3> key{n, m, l3};
  if (pnkb_undist_memo.count(key) == 0)
  {
    mpz_class binom1, binom2;
    mpz_bin_uiui(binom1.get_mpz_t(), n + 3 - l3, m + 2);
    mpz_bin_uiui(binom2.get_mpz_t(), n + 3, m + 3);
    mpq_class ret(binom1, binom2);
    pnkb_undist_memo.emplace(key, ret);
  }
  return pnkb_undist_memo[key];
}

MatrixCache& cached_matrices(const int n)
{
  if (cache.count(n) == 0)
  {
//    DEBUG1 << "moran eigensystem";
    const MoranEigensystem mei = compute_moran_eigensystem(n);
//    DEBUG1 << "moran eigensystem done";
    VectorXq D_subtend_above = VectorXq::LinSpaced(n, 1, n);
    D_subtend_above /= n + 1;

    VectorXq D_subtend_below = Eigen::Array<mpq_class, Eigen::Dynamic, 1>::Ones(n + 1) /
                               Eigen::Array<mpq_class, Eigen::Dynamic, 1>::LinSpaced(n + 1, 2, n + 2);
    D_subtend_below *= 2;

    MatrixXq Wnbj(n, n), P_dist(n + 1, n + 1), P_undist(n + 1, n);
    Wnbj.setZero();
    for (int b = 1; b < n + 1; ++b)
      for (int j = 2; j < n + 2; ++j)
        Wnbj(b - 1, j - 2) = calculate_Wnbj(n + 1, b, j);

    // P_dist(k, b) = probability of state (1, b) when there are k undistinguished lineages remaining
    P_dist.setZero();
    for (int k = 0; k < n + 1; ++k)
      for (int b = 1; b < n - k + 2; ++b)
        P_dist(k, b - 1) = pnkb_dist(n, k, b);

    // P_undist(k, b) = probability of state (0, b + 1) when there are k undistinguished lineages remaining
    P_undist.setZero();
    for (int k = 1; k < n + 1; ++k)
      for (int b = 1; b < n - k + 2; ++b)
        P_undist(k, b - 1) = pnkb_undist(n, k, b);

    VectorXq lsp = VectorXq::LinSpaced(n + 1, 2, n + 2);

//    DEBUG1 << "Eigen uses " << Eigen::nbThreads() << " for matmul";

    below_coeff bc = compute_below_coeffs(n);

    // This is too slow.
//    DEBUG1 << "X0";
    // parallel_matmul(Wnbj.transpose(),
    //         VectorXq::Ones(n) - D_subtend_above,
    //         mei.U.bottomRows(n),
    //         cache[n].X0);
    cache[n].X0 = (Wnbj.transpose() * (VectorXq::Ones(n) - D_subtend_above).asDiagonal() *
                   mei.U.bottomRows(n)).template cast<double>();
//    DEBUG1 << "X2";
    // parallel_matmul(Wnbj.transpose(),
    //         D_subtend_above,
    //         mei.U.reverse().topRows(n),
    //         cache[n].X2);
    cache[n].X2 = (Wnbj.transpose() * D_subtend_above.asDiagonal() *
                   mei.U.reverse().topRows(n)).template cast<double>();
//    DEBUG1 << "M0";
    // parallel_matmul(bc.coeffs,
    //         lsp.cwiseProduct((VectorXq::Ones(n + 1) - D_subtend_below)),
    //         P_undist,
    //         cache[n].M0);
    cache[n].M0 = (bc.coeffs * lsp.asDiagonal() * (VectorXq::Ones(n + 1) -
                                                   D_subtend_below).asDiagonal() * P_undist).template cast<double>();
//    DEBUG1 << "M1";
    cache[n].M1 = (bc.coeffs * lsp.asDiagonal() * D_subtend_below.asDiagonal() *
                   P_dist).template cast<double>();
    // parallel_matmul(bc.coeffs,
    //         lsp.cwiseProduct(D_subtend_below),
    //         P_dist,
    //         cache[n].M1);
    store_cache();
  }
  return cache.at(n);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// conditioned_sfs.cpp /////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
OnePopConditionedSFS<T>::OnePopConditionedSFS(int n) :
    n(n),
    mei(compute_moran_eigensystem(n)),
    mcache(cached_matrices(n)),
    Uinv_mp0(mei.Uinv.rightCols(n).template cast<double>()),
    Uinv_mp2(mei.Uinv.reverse().leftCols(n).template cast<double>())
{}

template <typename T>
std::vector<Matrix<T> > OnePopConditionedSFS<T>::compute_below(const PiecewiseConstantRateFunction<T> &eta) const
{
//  DEBUG1 << "compute below";
  const int M = eta.getHiddenStates().size() - 1;
  std::vector<Matrix<T> > csfs_below(M, Matrix<T>::Zero(3, n + 1));
  Matrix<T> tjj_below(M, n + 1);
  tjj_below.fill(eta.zero());
//  DEBUG1 << "tjj_double_integral below starts";
#pragma omp parallel for
  for (int m = 0; m < M; ++m)
    eta.tjj_double_integral_below(this->n, m, tjj_below);
//  DEBUG1 << "tjj_double_integral below finished";
//  DEBUG1 << "matrix products below (M0)";
  Matrix<T> M0_below = tjj_below * mcache.M0.template cast<T>();
//  DEBUG1 << "matrix products below (M1)";
  Matrix<T> M1_below = tjj_below * mcache.M1.template cast<T>();
//  DEBUG1 << "filling csfs_below";
  for (int m = 0; m < M; ++m)
  {
    csfs_below[m].fill(eta.zero());
    csfs_below[m].block(0, 1, 1, n) = M0_below.row(m);
    csfs_below[m].block(1, 0, 1, n + 1) = M1_below.row(m);
    CHECK_NAN(csfs_below[m]);
  }
//  DEBUG1 << "compute below finished";
  return csfs_below;
}

template <typename T>
std::vector<Matrix<T> > OnePopConditionedSFS<T>::compute_above(const PiecewiseConstantRateFunction<T> &eta) const
{
  const int M = eta.getHiddenStates().size() - 1;
  std::vector<Matrix<T> > C_above(M, Matrix<T>::Zero(n + 1, n));
  std::vector<Matrix<T> > csfs_above(M, Matrix<T>::Zero(3, n + 1));
//  DEBUG1 << "compute above";
#pragma omp parallel for
  for (int j = 2; j < n + 3; ++j)
    eta.tjj_double_integral_above(n, j, C_above);
  Matrix<T> tmp;

#pragma omp parallel for
  for (int m = 0; m < M; ++m)
  {
    csfs_above[m].fill(eta.zero());
    const Matrix<T> C0 = C_above[m].transpose();
    const Matrix<T> C2 = C_above[m].colwise().reverse().transpose();
    Vector<T> tmp0(this->mcache.X0.cols()), tmp2(this->mcache.X2.cols());
    tmp0.fill(eta.zero());
    for (int j = 0; j < this->mcache.X0.cols(); ++j)
    {
      std::vector<T> v;
      for (int i = 0; i < this->mcache.X0.rows(); ++i)
        v.push_back(this->mcache.X0(i, j) * C0(i, j));
      std::sort(v.begin(), v.end(), [] (T x, T y) { return std::abs(toDouble(x)) > std::abs(toDouble(y)); });
      tmp0(j) = doubly_compensated_summation(v);
    }
    csfs_above[m].block(0, 1, 1, n) = tmp0.transpose().lazyProduct(Uinv_mp0);
    tmp2.fill(eta.zero());
    for (int j = 0; j < this->mcache.X2.cols(); ++j)
    {
      std::vector<T> v;
      for (int i = 0; i < this->mcache.X2.rows(); ++i)
        v.push_back(this->mcache.X2(i, j) * C2(i, j));
      std::sort(v.begin(), v.end(), [] (T x, T y) { return std::abs(toDouble(x)) > std::abs(toDouble(y)); });
      tmp2(j) = doubly_compensated_summation(v);
    }
    csfs_above[m].block(2, 0, 1, n) = tmp2.transpose().lazyProduct(Uinv_mp2);
    CHECK_NAN(csfs_above[m]);
  }
  return csfs_above;
}

template <typename T>
std::vector<Matrix<T> > OnePopConditionedSFS<T>::compute(const PiecewiseConstantRateFunction<T> &eta)
{
//  DEBUG1 << "compute called";
  const int M = eta.getHiddenStates().size() - 1;
  std::vector<Matrix<T> > csfs_above = compute_above(eta);
  std::vector<Matrix<T> > csfs_below = compute_below(eta);
  std::vector<Matrix<T> > csfs(M, Matrix<T>::Zero(3, n + 1));
  for (int m = 0; m < M; ++m)
    csfs[m] = csfs_above[m] + csfs_below[m];
//  DEBUG1 << "compute finished";
  return csfs;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Matrix<adouble> sfs_cython(const int n, const ParameterVector p,
                           const double t1, const double t2, bool below_only)
{
  std::vector<double> hs{t1, t2};
  OnePopConditionedSFS<adouble> csfs(n);
  std::vector<Matrix<adouble> > v;
  PiecewiseConstantRateFunction<adouble> eta(p, hs);
  if (below_only)
    v = csfs.compute_below(eta);
  else
    v = csfs.compute(eta);
  return v.at(0);
}

/**
 * Converts the a and s vectors for vectors of type autodiff scalars, in a ParameterVector
 */
ParameterVector make_params(const std::vector<double>& a, const std::vector<double>& s) {
  ParameterVector ret;
  std::vector<adouble> r;
  for (const auto& aa : a) {
    r.emplace_back(aa);
  }
  ret.push_back(r);
  std::vector<adouble> cs;
  for (const auto& ss : s) {
    cs.emplace_back(ss);
  }
  ret.push_back(cs);
  return ret;
}

Matrix<double> raw_sfs(const std::vector<double>& a, const std::vector<double>& s, const int n, const double t1, const double t2, const bool below_only) {
  ParameterVector pv = make_params(a, s);
  Matrix<adouble> sfs = sfs_cython(n, pv, t1, t2, below_only);
  return sfs.unaryExpr([](adouble x){return x.value();});
}
