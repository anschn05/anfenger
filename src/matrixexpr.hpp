#ifndef FILE_MATEXPR
#define FILE_MATEXPR

#include <cassert>
#include <iostream>
#include <utility> // for std::declval

namespace ASC_bla
{
  template <typename T>
  class MatExpr
  { 
  public:
    auto derived() const { return static_cast<const T&> (*this); }
    size_t rows() const { return derived().rows(); }
    size_t cols() const { return derived().cols(); }
    auto operator() (size_t i, size_t j) const { return derived()(i,j); }
  };
  
  // ***************** Sum of two matrices *****************
  template <typename TA, typename TB>
  class SumMatExpr : public MatExpr<SumMatExpr<TA,TB>>
  {
    TA A;
    TB B;
  public:
    SumMatExpr (TA _A, TB _B) : A(_A), B(_B) { }
    auto operator() (size_t i, size_t j) const { return A(i,j) + B(i,j); }
    size_t rows() const { return A.rows(); }
    size_t cols() const { return A.cols(); }
  };

  template <typename TA, typename TB>
  auto operator+ (const MatExpr<TA> & A, const MatExpr<TB> & B)
  {
    assert (A.rows() == B.rows() && A.cols() == B.cols());
    return SumMatExpr(A.derived(), B.derived());
  }

  //  ***************** multiplication of two matrices (matrix product) *****************
  template <typename TA, typename TB>
  class MulMatExpr : public MatExpr<MulMatExpr<TA,TB>>
  {
    TA A;
    TB B;
  public:
    MulMatExpr (TA _A, TB _B) : A(_A), B(_B) { }

    // element (i,j) = sum_k A(i,k) * B(k,j)
    auto operator() (size_t i, size_t j) const {
      using ElemA = decltype(std::declval<TA>()(size_t{}, size_t{}));
      using ElemB = decltype(std::declval<TB>()(size_t{}, size_t{}));
      using SumT = decltype(std::declval<ElemA>() * std::declval<ElemB>());

      SumT sum = SumT();
      for (size_t k = 0; k < A.cols(); ++k)
        sum += A(i,k) * B(k,j);
      return sum;
    }

    size_t rows() const { return A.rows(); }
    size_t cols() const { return B.cols(); }
  };

  template <typename TA, typename TB>
  auto operator* (const MatExpr<TA> & A, const MatExpr<TB> & B)
  {
    assert(A.cols() == B.rows());
    return MulMatExpr(A.derived(), B.derived());
  }

   // ***************** Multiplication of vector and Matrix *****************
     template <typename TV, typename TM>
   class VecMatMulExpr : public MatExpr<VecMatMulExpr<TV,TM>>
   {
     TV v;
     TM m;
   public:
     VecMatMulExpr (TV _v, TM _m) : v(_v), m(_m) { }

     auto operator() (size_t i, size_t j) const {
       using ElemV = decltype(std::declval<TV>()(size_t{}));
       using ElemM = decltype(std::declval<TM>()(size_t{}, size_t{}));
       using ProdT = decltype(std::declval<ElemV>() * std::declval<ElemM>());

       ProdT sum = ProdT();
       for (size_t k = 0; k < m.rows(); ++k)
         sum += v(k) * m(k,j);
       return sum;
     }

     size_t rows() const { return 1; }
     size_t cols() const { return m.cols(); }
   };

  // output operator
  template <typename T>
  std::ostream & operator<< (std::ostream & ost, const MatExpr<T> & m)
  {
    for (size_t i = 0; i < m.rows(); ++i)
    {
      for (size_t j = 0; j < m.cols(); ++j)
      {
        ost << m(i,j);
        if (j + 1 < m.cols()) ost << ", ";
      }
      if (i + 1 < m.rows()) ost << "\n";
    }
    return ost;
  }

}


#endif
