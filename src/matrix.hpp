#ifndef FILE_MATRIX
#define FILE_MATRIX

#include <iostream>
#include <cassert>
#include "vector.hpp"

namespace ASC_bla
{
  template <typename T, typename TDIST = std::integral_constant<size_t,1>>
  class MatrixView
  {
  protected:
    T* m_data;
    size_t m_rows;
    size_t m_cols;
    TDIST m_rowdist; // Abstand zwischen Zeilen

  public:
    MatrixView() = default;
    MatrixView(const MatrixView&) = default;

    template<typename TDIST2>
    MatrixView(const MatrixView<T, TDIST2>& m2)
      : m_data(m2.data()), m_rows(m2.Rows()), m_cols(m2.Cols()), m_rowdist(m2.rowdist()) { }

    MatrixView(size_t rows, size_t cols, T* data)
      : m_data(data), m_rows(rows), m_cols(cols) { }

    MatrixView(size_t rows, size_t cols, TDIST rowdist, T* data)
      : m_data(data), m_rows(rows), m_cols(cols), m_rowdist(rowdist) { }

    // Zugriff
    T* data() const { return m_data; }
    size_t Rows() const { return m_rows; }
    size_t Cols() const { return m_cols; }
    auto rowdist() const { return m_rowdist; }

    // Elementzugriff (zeilenweise Speicherung)
    T& operator()(size_t i, size_t j) { return m_data[i*m_rowdist + j]; }
    const T& operator()(size_t i, size_t j) const { return m_data[i*m_rowdist + j]; }

    // Zeilen- und Spalten-Views
    auto Row(size_t i) const {
      assert(i < m_rows);
      return VectorView<T>(m_cols, m_data + i*m_rowdist);
    }

    auto Col(size_t j) const {
      assert(j < m_cols);
      return VectorView<T, size_t>(m_rows, m_rowdist, m_data + j);
    }

    // Teilmatrix (gleicher rowdist)
    auto Range(size_t first_row, size_t next_row, size_t first_col, size_t next_col) const {
      assert(first_row <= next_row && next_row <= m_rows);
      assert(first_col <= next_col && next_col <= m_cols);
      return MatrixView<T, TDIST>(
        next_row - first_row,
        next_col - first_col,
        m_rowdist,
        m_data + first_row*m_rowdist + first_col
      );
    }

    MatrixView& operator=(T val) {
      for (size_t i = 0; i < m_rows; i++)
        for (size_t j = 0; j < m_cols; j++)
          (*this)(i,j) = val;
      return *this;
    }

    template<typename TB>
    MatrixView& operator=(const TB& other) {
      assert(m_rows == other.Rows() && m_cols == other.Cols());
      for (size_t i = 0; i < m_rows; i++)
        for (size_t j = 0; j < m_cols; j++)
          (*this)(i,j) = other(i,j);
      return *this;
    }
  };

/*
  template <typename T>
  class Matrix : public MatrixView<T>
  {
    using BASE = MatrixView<T>;
    using BASE::m_data;
    using BASE::m_rows;
    using BASE::m_cols;
    using BASE::m_rowdist;

  public:
    Matrix(size_t rows, size_t cols)
      : MatrixView<T>(rows, cols, cols, new T[rows * cols]){ }

    Matrix(const Matrix& m)
      : Matrix(m.Rows(), m.Cols())
    {
      *this = m;
    }

    Matrix(Matrix&& m)
      : MatrixView<T>(0, 0, 0, nullptr)
    {
      std::swap(m_data, m.m_data);
      std::swap(m_rows, m.m_rows);
      std::swap(m_cols, m.m_cols);
      std::swap(m_rowdist, m.m_rowdist);
    }

    template<typename TB>
    Matrix(const TB& m)
      : Matrix(m.Rows(), m.Cols())
    {
      *this = m;
    }

    ~Matrix() { delete[] m_data; }

    using BASE::operator=;

    Matrix& operator=(const Matrix& m2) {
      assert(m_rows == m2.m_rows && m_cols == m2.m_cols);
      for (size_t i = 0; i < m_rows; i++)
        for (size_t j = 0; j < m_cols; j++)
          (*this)(i,j) = m2(i,j);
      return *this;
    }

    Matrix& operator=(Matrix&& m2) {
      std::swap(m_data, m2.m_data);
      std::swap(m_rows, m2.m_rows);
      std::swap(m_cols, m2.m_cols);
      std::swap(m_rowdist, m2.m_rowdist);
      return *this;
    }
  };


  // Ausgabeoperator
  template <typename T>
  std::ostream& operator<<(std::ostream& ost, const Matrix<T>& M)
  {
    for (size_t i = 0; i < M.Rows(); i++) {
      for (size_t j = 0; j < M.Cols(); j++) {
        ost << M(i,j);
        if (j + 1 < M.Cols()) ost << ", ";
      }
      if (i + 1 < M.Rows()) ost << "\n";
    }
    return ost;
  }
}

#endif
*/

/*
#ifndef FILE_MATRIX
#define FILE_MATRIX

#include <iostream>
#include <utility>   // für std::swap
#include <cassert>   // für assert

namespace ASC_bla
{
*/
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


