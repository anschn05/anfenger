#ifndef FILE_MATRIX
#define FILE_MATRIX

#include <iostream>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include "matrixexpr.hpp"
#include "vector.hpp"

namespace ASC_bla
{

// =======================================
// MatrixView – non-owning view + Expr interface
// =======================================
template <typename T, typename TDIST = size_t>
class MatrixView : public MatExpr<MatrixView<T, TDIST>>
{
protected:
    T* m_data;
    size_t m_rows;
    size_t m_cols;
    TDIST m_rowdist;

public:
    MatrixView() = default;
    MatrixView(const MatrixView&) = default;

    template<typename TDIST2>
    MatrixView(const MatrixView<T, TDIST2>& m2)
        : m_data(m2.data()), m_rows(m2.Rows()), m_cols(m2.Cols()), m_rowdist(m2.rowdist()) { }

    MatrixView(size_t rows, size_t cols, T* data)
        : m_data(data), m_rows(rows), m_cols(cols), m_rowdist(cols) { }

    MatrixView(size_t rows, size_t cols, TDIST rowdist, T* data)
        : m_data(data), m_rows(rows), m_cols(cols), m_rowdist(rowdist) { }

    // Accessors
    T* data() const { return m_data; }
    size_t rows() const { return m_rows; }
    size_t cols() const { return m_cols; }
    size_t Rows() const { return m_rows; }
    size_t Cols() const { return m_cols; }
    auto rowdist() const { return m_rowdist; }

    // Element access
    T& operator()(size_t i, size_t j) { return m_data[i*m_rowdist + j]; }
    const T& operator()(size_t i, size_t j) const { return m_data[i*m_rowdist + j]; }

    // Row/Column views
    auto Row(size_t i) const {
        assert(i < m_rows);
        return VectorView<T>(m_cols, m_data + i*m_rowdist);
    }

    auto Col(size_t j) const {
        assert(j < m_cols);
        return VectorView<T, size_t>(m_rows, m_rowdist, m_data + j);
    }

    // Submatrix range
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

    // Assign scalar
    MatrixView& operator=(T val) {
        for (size_t i = 0; i < m_rows; i++)
            for (size_t j = 0; j < m_cols; j++)
                (*this)(i,j) = val;
        return *this;
    }

    // Assign from another matrix/expression
    template<typename TB>
    MatrixView& operator=(const TB& other) {
        assert(m_rows == other.Rows() && m_cols == other.Cols());
        for (size_t i = 0; i < m_rows; i++)
            for (size_t j = 0; j < m_cols; j++)
                (*this)(i,j) = other(i,j);
        return *this;
    }
};


// =======================================
// Matrix – owning version
// =======================================
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
        : MatrixView<T>(rows, cols, cols, new T[rows * cols]) { }

    Matrix(const Matrix& m)
        : Matrix(m.Rows(), m.Cols()) { *this = m; }

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


// =======================================
// Ausgabeoperator
// =======================================
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


// =======================================
// Numerisch stabile Pivot-Suche
// =======================================
template <typename T>
size_t FindPivot(const Matrix<T>& mat, size_t col, size_t start_row, double tolerance = 1e-12)
{
    size_t pivot_row = start_row;
    T max_val = std::abs(mat(start_row, col));

    for (size_t i = start_row + 1; i < mat.Rows(); ++i) {
        T current_val = std::abs(mat(i, col));
        if (current_val > max_val) {
            max_val = current_val;
            pivot_row = i;
        }
    }

    if (max_val < tolerance)
        throw std::runtime_error("Matrix is singular - zero column detected");

    return pivot_row;
}


// =======================================
// Matrix-Inversion mit Pivoting (Gauß–Jordan)
// =======================================
template <typename T>
Matrix<T> Inverse(const Matrix<T>& m)
{
    assert(m.Rows() == m.Cols() && "Matrix must be square to compute its inverse");
    size_t n = m.Rows();
    Matrix<T> augmented(n, 2 * n);

    // augmented = [A | I]
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) {
            augmented(i, j) = m(i, j);
            augmented(i, j + n) = (i == j) ? T(1) : T(0);
        }

    // Gaussian elimination with partial pivoting
    for (size_t i = 0; i < n; ++i)
    {
        size_t pivot_row = i;
        if (std::abs(augmented(i,i)) < 1e-12)
            pivot_row = FindPivot(augmented, i, i);

        if (pivot_row != i)
            for (size_t j = 0; j < 2 * n; ++j)
                std::swap(augmented(i,j), augmented(pivot_row,j));

        T pivot = augmented(i, i);
        assert(std::abs(pivot) >= 1e-12 && "Matrix is singular and cannot be inverted");

        for (size_t j = 0; j < 2 * n; ++j)
            augmented(i, j) /= pivot;

        for (size_t k = 0; k < n; ++k)
            if (k != i) {
                T factor = augmented(k, i);
                for (size_t j = 0; j < 2 * n; ++j)
                    augmented(k, j) -= factor * augmented(i, j);
            }
    }

    // Extract right block [I | A⁻¹]
    Matrix<T> inverse(n, n);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            inverse(i, j) = augmented(i, j + n);

    return inverse;
}

} // namespace ASC_bla

#endif
