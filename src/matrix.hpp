
/*#ifndef FILE_MATRIX
#define FILE_MATRIX

#include <iostream>
#include <utility>   // für std::swap
#include <cassert>   // für assert

namespace ASC_bla
{

  template <typename T>
  class Matrix
  {
    size_t nrows;
    size_t ncols;
    T* data; // ?

  public:
    // Konstruktor
    Matrix(size_t _nrows, size_t _ncols)
      : nrows(_nrows), ncols(_ncols), data(new T[_nrows * _ncols]) { }  // ?

    // Copy-Konstruktor
    Matrix(const Matrix& m)
      : Matrix(m.nrows, m.ncols)
    {
      *this = m;
    }

    // Move-Konstruktor
    Matrix(Matrix&& m)
      : nrows(0), ncols(0), data(nullptr)
    {
      std::swap(nrows, m.nrows);
      std::swap(ncols, m.ncols);
      std::swap(data, m.data);
    }

    // Destruktor
    ~Matrix() { delete[] data; }

    // Kopierzuweisung
    Matrix& operator=(const Matrix& m2)
    {
      assert(nrows == m2.nrows && ncols == m2.ncols && "Matrix dimensions must match");
      for (size_t i = 0; i < nrows * ncols; i++)
        data[i] = m2.data[i];
      return *this;
    }

    // Move-Zuweisung
    Matrix& operator=(Matrix&& m2)
    {
      std::swap(nrows, m2.nrows);
      std::swap(ncols, m2.ncols);
      std::swap(data, m2.data);
      return *this;
    }

    // Zugriff mit (i,j)
    T& operator()(size_t i, size_t j) { return data[i * ncols + j]; }
    const T& operator()(size_t i, size_t j) const { return data[i * ncols + j]; }

    // Zugriff linear (wie Vector)
    T& operator()(size_t idx) { return data[idx]; }
    const T& operator()(size_t idx) const { return data[idx]; }

    // Getter
    size_t Rows() const { return nrows; }
    size_t Cols() const { return ncols; }
    size_t Size() const { return nrows * ncols; }
  };

  // Addition
  template <typename T>
  Matrix<T> operator+(const Matrix<T>& a, const Matrix<T>& b)
  {
    assert(a.Rows() == b.Rows() && a.Cols() == b.Cols() && "Matrix dimensions must match");
    Matrix<T> sum(a.Rows(), a.Cols());
    for (size_t i = 0; i < a.Size(); i++)
      sum(i) = a(i) + b(i);
    return sum;
  }

  // Subtraktion
  template <typename T>
  Matrix<T> operator-(const Matrix<T>& a, const Matrix<T>& b)
  {
    assert(a.Rows() == b.Rows() && a.Cols() == b.Cols() && "Matrix dimensions must match");
    Matrix<T> diff(a.Rows(), a.Cols());
    for (size_t i = 0; i < a.Size(); i++)
      diff(i) = a(i) - b(i);
    return diff;
  }

  // Ausgabeoperator
  template <typename T>
  std::ostream& operator<<(std::ostream& ost, const Matrix<T>& m)
  {
    for (size_t i = 0; i < m.Rows(); i++)
    {
      for (size_t j = 0; j < m.Cols(); j++)
      {
        ost << m(i, j);
        if (j + 1 < m.Cols()) ost << ", ";
      }
      if (i + 1 < m.Rows()) ost << "\n";
    }
    return ost;
  }

    // Matrix-Multiplikation
    template <typename T>
    Matrix<T> operator*(const Matrix<T>& a, const Matrix<T>& b)
    {
      assert(a.Cols() == b.Rows() && "Matrix dimensions must be compatible for multiplication");
      Matrix<T> result(a.Rows(), b.Cols());
      for (size_t i = 0; i < a.Rows(); ++i)
        for (size_t j = 0; j < b.Cols(); ++j)
        {
          T sum = T();
          for (size_t k = 0; k < a.Cols(); ++k)
            sum += a(i, k) * b(k, j);
          result(i, j) = sum;
        }
      return result;
    }

    //Matrix Inverse
    template <typename T>
    Matrix<T> Inverse(const Matrix<T>& m)
      {
        assert(m.Rows() == m.Cols() && "Matrix must be square to compute its inverse");
        size_t n = m.Rows();
        Matrix<T> augmented(n, 2 * n);

        // New Matrix to use Gauß-Algorithm, with id
        for (size_t i = 0; i < n; ++i)
          for (size_t j = 0; j < n; ++j)
          {
            augmented(i, j) = m(i, j);
            augmented(i, j + n) = (i == j) ? T(1) : T(0);
          }

        // Gauß-Algorithmus
        for (size_t i = 0; i < n; ++i)
        {
          // Look for pivot (non zero in current column)
          T pivot = augmented(i, i);
          assert(pivot != T(0) && "Matrix is singular and cannot be inverted");

          for (size_t j = 0; j < 2 * n; ++j)
            augmented(i, j) /= pivot;

          // elim other entries in curr col
          for (size_t k = 0; k < n; ++k)
          {
            if (k != i)
            {
              T factor = augmented(k, i);
              for (size_t j = 0; j < 2 * n; ++j)
                augmented(k, j) -= factor * augmented(i, j);
            }
          }
        }

        // end gauss by extracting inv
        Matrix<T> inverse(n, n);
        for (size_t i = 0; i < n; ++i)
          for (size_t j = 0; j < n; ++j)
            inverse(i, j) = augmented(i, j + n);

        return inverse;
      }

} // namespace ASC_bla

#endif
*/

#ifndef FILE_MATRIX
#define FILE_MATRIX

#include <cassert>
#include <cstddef>
#include <utility>
#include <iostream>

#include "vecexpr.hpp"   // wie bei dir
#include "vector.hpp"    // dein Header mit Vector/VectorView (oben gepostet)

namespace ASC_bla
{

  // MatrixView: nicht-ownend, row-major; col-dist = 1 implizit.
  // Optionaler rowdist (Leading Dimension) für Subviews/Slices.
  template <typename T, typename TROWDIST = std::integral_constant<size_t, 0>>
  class MatrixView /* duck-typed MatExpr */ {
  protected:
    T*      m_data = nullptr;
    size_t  m_rows = 0;
    size_t  m_cols = 0;
    // Wenn TROWDIST == integral_constant<...,0>, wird rowdist zur Laufzeit = m_cols gesetzt.
    TROWDIST m_rowdist{};

    size_t effective_rowdist() const {
      if constexpr (std::is_integral_v<TROWDIST>)
        return static_cast<size_t>(m_rowdist);
      else
        return m_rowdist; // falls jemand size_t übergibt
    }
    size_t rowdist_runtime() const {
      if constexpr (std::is_same_v<TROWDIST, std::integral_constant<size_t,0>>)
        return m_cols;            // contiguous default
      else
        return effective_rowdist();
    }

  public:
    MatrixView() = default;
    MatrixView(const MatrixView&) = default;

    template <typename TROWDIST2>
    MatrixView(const MatrixView<T, TROWDIST2>& mv)
      : m_data(mv.data()), m_rows(mv.rows()), m_cols(mv.cols()), m_rowdist(mv.rowdist()) {}

    // Contiguous view (rowdist = cols)
    MatrixView(size_t r, size_t c, T* data)
      : m_data(data), m_rows(r), m_cols(c) {}

    // View mit explizitem rowdist (ld)
    MatrixView(size_t r, size_t c, TROWDIST rowdist, T* data)
      : m_data(data), m_rows(r), m_cols(c), m_rowdist(rowdist) {}

    // Zuweisung aus "MatExpr"-ähnlichem Typ (duck typing: .rows(), .cols(), (i,j))
    template <typename MB>
    MatrixView& operator=(const MB& rhs) {
      assert(m_rows == rhs.rows() && m_cols == rhs.cols());
      const size_t ld = rowdist_runtime();
      for (size_t i = 0; i < m_rows; ++i)
        for (size_t j = 0; j < m_cols; ++j)
          m_data[i*ld + j] = rhs(i, j);
      return *this;
    }

    // Skalarfüllung
    MatrixView& operator=(const T& s) {
      const size_t ld = rowdist_runtime();
      for (size_t i = 0; i < m_rows; ++i)
        for (size_t j = 0; j < m_cols; ++j)
          m_data[i*ld + j] = s;
      return *this;
    }

    // Accessors
    T*      data() const { return m_data; }
    size_t  rows() const { return m_rows; }
    size_t  cols() const { return m_cols; }
    size_t  rowdist() const { return rowdist_runtime(); }

    // Elementzugriff
    T& operator()(size_t i, size_t j) { return m_data[i*rowdist_runtime() + j]; }
    const T& operator()(size_t i, size_t j) const { return m_data[i*rowdist_runtime() + j]; }

    // Row-/Col-Views (nutzen deinen VectorView)
    auto row(size_t i) const {
      assert(i < m_rows);
      return VectorView<T, size_t>(m_cols, size_t(1), m_data + i*rowdist_runtime());
    }
    auto col(size_t j) const {
      assert(j < m_cols);
      return VectorView<T, size_t>(m_rows, rowdist_runtime(), m_data + j);
    }

    // Block-View (Subview)
    auto block(size_t i0, size_t j0, size_t r, size_t c) const {
      assert(i0 + r <= m_rows && j0 + c <= m_cols);
      return MatrixView<T, size_t>(r, c, rowdist_runtime(), m_data + i0*rowdist_runtime() + j0);
    }
  };


  // Owning Matrix (wie dein Vector)
  template <typename T>
  class Matrix : public MatrixView<T> {
    using BASE = MatrixView<T>;
    using BASE::m_data;  using BASE::m_rows;  using BASE::m_cols;

  public:
    Matrix() = default;

    Matrix(size_t r, size_t c)
      : MatrixView<T>(r, c, new T[r*c]) {}

    Matrix(const Matrix& other)
      : Matrix(other.rows(), other.cols()) {
      *this = other;
    }

    Matrix(Matrix&& other) noexcept
      : MatrixView<T>(0, 0, nullptr) {
      std::swap(m_data, other.m_data);
      std::swap(m_rows, other.m_rows);
      std::swap(m_cols, other.m_cols);
    }

    template <typename MB>
    Matrix(const MB& rhs)
      : Matrix(rhs.rows(), rhs.cols()) {
      *this = rhs;
    }

    ~Matrix() { delete [] m_data; }

    using BASE::operator=; // erbt Zuweisungen (Skalar + MatExpr)

    Matrix& operator=(const Matrix& rhs) {
      assert(m_rows == rhs.rows() && m_cols == rhs.cols());
      const size_t ld = this->rowdist();
      for (size_t i = 0; i < m_rows; ++i)
        for (size_t j = 0; j < m_cols; ++j)
          m_data[i*ld + j] = rhs(i, j);
      return *this;
    }

    Matrix& operator=(Matrix&& rhs) noexcept {
      std::swap(m_data, rhs.m_data);
      std::swap(m_rows, rhs.m_rows);
      std::swap(m_cols, rhs.m_cols);
      return *this;
    }
  };


  // Einfache Operationen (ohne temporäre Ausdrucks-Templates, straight loops)
  template <typename T, typename A, typename B>
  struct MatAdd {
    const A& a; const B& b;
    size_t rows() const { return a.rows(); }
    size_t cols() const { return a.cols(); }
    T operator()(size_t i, size_t j) const { return a(i,j) + b(i,j); }
  };

  template <typename T>
  Matrix<T> operator+(const MatrixView<T>& a, const MatrixView<T>& b) {
    assert(a.rows()==b.rows() && a.cols()==b.cols());
    Matrix<T> out(a.rows(), a.cols());
    out = MatAdd<T, MatrixView<T>, MatrixView<T>>{a, b};
    return out;
  }

  template <typename T>
  Matrix<T> operator-(const MatrixView<T>& a, const MatrixView<T>& b) {
    assert(a.rows()==b.rows() && a.cols()==b.cols());
    Matrix<T> out(a.rows(), a.cols());
    struct M { const MatrixView<T>& a; const MatrixView<T>& b;
      size_t rows() const { return a.rows(); }
      size_t cols() const { return a.cols(); }
      T operator()(size_t i, size_t j) const { return a(i,j) - b(i,j); }
    } minus{a,b};
    out = minus;
    return out;
  }

  // Matrix * Vector -> Vector
  template <typename T>
  Vector<T> operator*(const MatrixView<T>& A, const VectorView<T>& x) {
    assert(A.cols() == x.size());
    Vector<T> y(A.rows());
    const size_t ld = A.rowdist();
    for (size_t i = 0; i < A.rows(); ++i) {
      T acc{}; // 0-initialisiert
      for (size_t j = 0; j < A.cols(); ++j)
        acc += A.data()[i*ld + j] * x(j);
      y(i) = acc;
    }
    return y;
  }

  // Ausgabe
  template <typename T>
  std::ostream& operator<<(std::ostream& os, const MatrixView<T>& M) {
    const size_t ld = M.rowdist();
    for (size_t i = 0; i < M.rows(); ++i) {
      if (M.cols() > 0) os << M.data()[i*ld + 0];
      for (size_t j = 1; j < M.cols(); ++j)
        os << ", " << M.data()[i*ld + j];
      if (i+1 < M.rows()) os << '\n';
    }
    return os;
  }

} // namespace ASC_bla

#endif // FILE_MATRIX
