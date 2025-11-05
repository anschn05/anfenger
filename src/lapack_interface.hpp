#ifndef FILE_LAPACK_INTERFACE_HPP
#define FILE_LAPACK_INTERFACE_HPP

#include <iostream>
#include <string>

#include "vector.hpp"
#include "matrix.hpp"



#include <complex>
#include <vector>
#include <algorithm>
#include <stdexcept>

typedef int integer;
typedef integer logical;
typedef float real;
typedef double doublereal;
typedef std::complex<float> singlecomplex;
typedef std::complex<double> doublecomplex;

// Windows SDK defines VOID in the file WinNT.h
#ifndef VOID
typedef void VOID;
#endif
typedef int ftnlen;
typedef int L_fp;  


extern "C" {
#include <clapack.h>
}


namespace ASC_bla
{


  // BLAS-1 functions:

  /*
    int daxpy_(integer *n, doublereal *da, doublereal *dx, 
    integer *incx, doublereal *dy, integer *incy);
  */
  // y += alpha x
  template <typename SX, typename SY>
  void addVectorLapack (double alpha, VectorView<double,SX> x, VectorView<double,SY> y)
  {
    integer n = x.size();
    integer incx = x.dist();
    integer incy = y.dist();
    daxpy_ (&n, &alpha, &x(0),  &incx, &y(0), &incy);
  }
  
  
  // BLAS-2 functions:

  // BLAS-3 functions:
  
  // BLAS-3 functions (dgemm) - wrapper that adapts our MatrixView layout
  // We rely on the declaration provided by the included clapack.h to declare dgemm_.
  // Some LAPACK/BLAS headers declare dgemm_ with slightly different linkage/signatures;
  // if you get warnings about conflicting declarations, prefer using the system-provided
  // header (clapack.h) only and remove local forward-declarations.

  // c = a*b using LAPACK/BLAS dgemm_. This version works with our MatrixView
  // which encodes storage via rowdist(): if rowdist()==1 it's column-major,
  // otherwise it's row-major (contiguous row-major usually rowdist==cols()).
  template <typename TA, typename TB, typename TC>
  void multMatMatLapack (const MatrixView<double, TA> & a,
                         const MatrixView<double, TB> & b,
                         MatrixView<double, TC> & c)
  {
    // If rowdist == 1 => column-major storage -> no transpose for dgemm
    char transa = (a.rowdist() == 1) ? 'N' : 'T';
    char transb = (b.rowdist() == 1) ? 'N' : 'T';

    integer m = static_cast<integer>(c.rows());
    integer n = static_cast<integer>(c.cols());
    integer k = static_cast<integer>( (transa == 'N') ? a.cols() : a.rows() );

    double alpha = 1.0;
    double beta = 0.0;
    integer lda = static_cast<integer>(std::max<size_t>(a.rowdist(), 1));
    integer ldb = static_cast<integer>(std::max<size_t>(b.rowdist(), 1));
    integer ldc = static_cast<integer>(std::max<size_t>(c.rowdist(), 1));

    // call Fortran dgemm (note Fortran expects char pointers)
    dgemm_(&transa, &transb, &m, &n, &k, &alpha,
           a.data(), &lda, b.data(), &ldb, &beta, c.data(), &ldc);
  }
  

  

  // LAPACK-backed LU factorization wrapper.
  // This class stores the LU factorization in a column-major buffer (as LAPACK expects)
  // while the project's Matrix/MatrixView types are row-major. We copy into a
  // column-major array for LAPACK calls and convert back when needed.
  class LapackLU {
    std::vector<doublereal> a_col; // column-major storage of LU factors
    std::vector<integer> ipiv;
    integer n = 0;
  public:
    // Factorize a (copy of input Matrix<double>)
    LapackLU(const Matrix<double> & A)
    {
      n = static_cast<integer>(A.rows());
      if (n <= 0) return;
      a_col.assign(n * n, 0.0);
      ipiv.assign(n, 0);

      // copy A (row-major) into column-major buffer a_col
      for (integer i = 0; i < n; ++i)
        for (integer j = 0; j < n; ++j)
          a_col[j*n + i] = A(static_cast<size_t>(i), static_cast<size_t>(j));

      integer info = 0;
      // dgetrf_(m,n,a,lda,ipiv,info)
      dgetrf_(&n, &n, a_col.data(), &n, ipiv.data(), &info);
      if (info < 0)
        throw std::runtime_error("dgetrf_: illegal argument");
      if (info > 0)
        throw std::runtime_error("dgetrf_: matrix is singular; U(" + std::to_string(info) + ") is zero");
    }

    // Solve A x = b with the computed LU; b is overwritten with the solution
  void solve(VectorView<double> b) {
      if (n == 0) return;
      if (static_cast<integer>(b.size()) != n)
        throw std::runtime_error("LapackLU::solve: size mismatch");

      char trans = 'N';
      integer nrhs = 1;
      integer ldb = static_cast<integer>(b.size());
      integer info = 0;
      // dgetrs_(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
      // Note: dgetrs_ expects non-const a, but does not modify LU; cast away const is OK here
      dgetrs_(&trans, &n, &nrhs, const_cast<doublereal*>(a_col.data()), &n, const_cast<integer*>(ipiv.data()), b.data(), &ldb, &info);
      if (info != 0)
        throw std::runtime_error("dgetrs_ failed with info=" + std::to_string(info));
    }

    // Compute inverse of the original matrix (returns a Matrix<double> in row-major order)
  Matrix<double> inverse() {
      if (n == 0) return Matrix<double>();
      // make a copy since dgetri overwrites the matrix
      std::vector<doublereal> a_inv = a_col;
      integer info = 0;
      // query optimal worksize
      doublereal work_query = 0.0;
      integer lwork = -1;
      dgetri_(&n, a_inv.data(), &n, const_cast<integer*>(ipiv.data()), &work_query, &lwork, &info);
      if (info != 0)
        throw std::runtime_error("dgetri_ workspace query failed");
      lwork = static_cast<integer>(work_query);
      std::vector<doublereal> work(std::max(1, static_cast<int>(lwork)));
      dgetri_(&n, a_inv.data(), &n, const_cast<integer*>(ipiv.data()), work.data(), &lwork, &info);
      if (info != 0)
        throw std::runtime_error("dgetri_ failed with info=" + std::to_string(info));

      // a_inv now holds inverse in column-major; convert to Matrix<double> row-major
      Matrix<double> M(static_cast<size_t>(n), static_cast<size_t>(n));
      for (integer i = 0; i < n; ++i)
        for (integer j = 0; j < n; ++j)
          M(static_cast<size_t>(i), static_cast<size_t>(j)) = a_inv[j*n + i];
      return M;
    }
  };
  
}


#endif
