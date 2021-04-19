// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef PREPAREDECODING_SMCPP_HPP
#define PREPAREDECODING_SMCPP_HPP

#include <iostream>
#include <fstream>
#include <random>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/AutoDiff>
#include "mpreal.h"

#include <fmt/format.h>
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ranges.h>
#include <fmt/ostream.h>

#ifdef NO_CHECK_NAN
#define CHECK_NAN(x)
#define CHECK_NAN_OR_NEGATIVE(x)
#else
#define CHECK_NAN(x) check_nan(x, __FILE__, __LINE__)
#define CHECK_NAN_OR_NEGATIVE(x) check_nan(x, __FILE__, __LINE__); check_negative(x, __FILE__, __LINE__);
#endif

#define ERROR std::cout << "smcpp: "
#define DEBUG1 std::cout << "smcpp: "

template <typename T> using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <typename T> using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <typename T, size_t P> using FixedVector = Eigen::Matrix<T, P, 1>;

typedef double adouble_base_type;
typedef Eigen::Matrix<adouble_base_type, Eigen::Dynamic, 1> adouble_t;
typedef Eigen::AutoDiffScalar<adouble_t> adouble;
typedef std::vector<std::vector<adouble>> ParameterVector;

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
void smcpp_init_cache();


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

Matrix<adouble> sfs_cython(const int, const ParameterVector, const double, const double, bool);

class PdModel {
private:
  std::vector<double> m_a;
  std::vector<double> m_s;

public:
  PdModel(std::vector<double> a, std::vector<double> s, const double N0) : m_a(std::move(a)), m_s(std::move(s)) {
  }

  [[nodiscard]] const std::vector<double>& getA() const {
    return m_a;
  }
  [[nodiscard]] const std::vector<double>& getS() const {
    return m_s;
  }
};

ParameterVector make_params(const std::vector<double>& a, const std::vector<double>& s);
ParameterVector make_params_from_model(const PdModel& model);
Matrix<double> raw_sfs(const PdModel& model, int n, double t1, double t2, bool below_only=false);


#endif // PREPAREDECODING_SMCPP_HPP
