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

} // namespace ASC_bla

#endif