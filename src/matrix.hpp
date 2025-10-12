#ifndef FILE_MATRIX
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