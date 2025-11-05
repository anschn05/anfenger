#include <iostream>
#include <chrono>
#include <cmath>

#include "matrix.hpp"
#include "lapack_interface.hpp"

using namespace ASC_bla;
using namespace std::chrono;

int main()
{
  const size_t N = 200; // Testgröße (kleiner zum Debug, größer für Messungen)
  Matrix<double> A(N,N), B(N,N);
  // Fill matrices (example: deterministic values)
  for (size_t i=0;i<N;i++)
    for (size_t j=0;j<N;j++) {
      A(i,j) = double(i + 1) * 0.001 + double(j+1)*0.002;
      B(i,j) = double(i + 1) * 0.003 - double(j+1)*0.001;
    }

  Matrix<double> C_naive(N,N), C_blas(N,N);

  // naive multiplication
  auto t0 = high_resolution_clock::now();
  for (size_t i=0;i<N;i++)
    for (size_t j=0;j<N;j++) {
      double s = 0.0;
      for (size_t k=0;k<N;k++)
        s += A(i,k) * B(k,j);
      C_naive(i,j) = s;
    }
  auto t1 = high_resolution_clock::now();
  double naive_time = duration_cast<duration<double>>(t1-t0).count();

  // LAPACK/BLAS multiplication via wrapper
  t0 = high_resolution_clock::now();
  // Note: our wrapper expects MatrixView arguments; Matrix<T> derives MatrixView<T>, so this deduces template args
  ASC_bla::multMatMatLapack(A, B, C_blas);
  t1 = high_resolution_clock::now();
  double lapack_time = duration_cast<duration<double>>(t1-t0).count();

  // compare results
  double maxdiff = 0.0;
  for (size_t i=0;i<N;i++)
    for (size_t j=0;j<N;j++) {
      double d = std::abs(C_naive(i,j) - C_blas(i,j));
      if (d > maxdiff) maxdiff = d;
    }

  std::cout << "N = " << N << "\n";
  std::cout << "naive time:  " << naive_time << " s\n";
  std::cout << "lapack time: " << lapack_time << " s\n";
  std::cout << "max abs difference: " << maxdiff << "\n";

  // basic check with tolerance
  double tol = 1e-9;
  if (maxdiff < tol) std::cout << "OK: results match within tol=" << tol << "\n";
  else std::cout << "WARNING: difference exceeds tol\n";

  return 0;
}