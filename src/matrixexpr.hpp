#ifndef FILE_MATEXPR
#define FILE_MATEXPR

#include <cassert>
#include <iostream>

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

  //  ***************** multiplication of two matrices *****************
  template <typename TA, typename TB>
    class MulMatExpr : public MatExpr<MulMatExpr<TA,TB>>
  {
    TA A;
    TB B;
  public:
  //MulMatExpr (TA _A, TB _B) : A(_A), B(_B) { }
   // auto operator() (size_t i, size_t j) const { return A(i,j) * B(i,j); }
   // size_t rows() const { return A.rows(); }
   // size_t cols() const { return A.cols(); }
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
