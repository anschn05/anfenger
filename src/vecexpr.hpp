#ifndef FILE_VECEXPR
#define FILE_VECEXPR

#include <cassert>
#include <utility>

namespace ASC_bla
{
  template <typename T>
  class VecExpr
  {
  public:
    auto derived() const { return static_cast<const T&>(*this); }
    size_t size() const { return derived().size(); }
    auto operator() (size_t i) const { return derived()(i); }
  };

  // Sum of two vectors
  template <typename TA, typename TB>
  class SumVecExpr : public VecExpr<SumVecExpr<TA,TB>>
  {
    TA A; TB B;
  public:
    SumVecExpr (TA _A, TB _B) : A(_A), B(_B) { }
    auto operator() (size_t i) const { return A(i) + B(i); }
    size_t size() const { return A.size(); }
  };

  template <typename TA, typename TB>
  auto operator+ (const VecExpr<TA> & A, const VecExpr<TB> & B)
  {
    assert (A.size() == B.size());
  return SumVecExpr<TA,TB>(A.derived(), B.derived());
  }

  // scalar * vector
  template <typename TV>
  class ScalMulVecExpr : public VecExpr<ScalMulVecExpr<TV>>
  {
    double s;
    TV v;
  public:
    ScalMulVecExpr(double _s, TV _v) : s(_s), v(_v) { }
    auto operator() (size_t i) const { return s * v(i); }
    size_t size() const { return v.size(); }
  };

  template <typename TV>
  auto operator* (double s, const VecExpr<TV> & v)
  {
  return ScalMulVecExpr<TV>(s, v.derived());
  }

  template <typename TV>
  auto operator* (const VecExpr<TV> & v, double s)
  {
  return ScalMulVecExpr<TV>(s, v.derived());
  }

  // dot product
  template <typename TA, typename TB>
  auto dot(const VecExpr<TA> & A, const VecExpr<TB> & B)
  {
    assert(A.size() == B.size());
    using ElemA = decltype(std::declval<TA>()(size_t{}));
    using ElemB = decltype(std::declval<TB>()(size_t{}));
    using SumT = decltype(std::declval<ElemA>() * std::declval<ElemB>());
    SumT sum = SumT();
    for (size_t i = 0; i < A.size(); ++i)
      sum += A(i) * B(i);
    return sum;
  }

}

#endif
