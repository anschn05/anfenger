#include <iostream>
#include <cmath>

#include "../src/matrix.hpp"
#include "../src/vector.hpp"
#include "../src/lapack_interface.hpp"

using namespace ASC_bla;
using std::cout; using std::endl;

int main()
{
  const size_t N = 5; // matrix size (>=5)

  // Build a non-singular test matrix
  Matrix<double> A(N,N);
  // make A diagonally dominant to avoid singular matrices
  for (size_t i=0;i<N;i++)
    for (size_t j=0;j<N;j++)
      A(i,j) = (i==j) ? (10.0 + double(i)) : 1.0;

  // Right-hand side b (set so solution exists)
  Vector<double> b(N);
  for (size_t i=0;i<N;i++) b(i) = double(i) + 1.0;

  cout << "A =\n" << A << "\n";
  cout << "b = " << b << "\n";

  // Factorize A using LAPACK
  LapackLU lu(A);

  // Solve system A x = b
  Vector<double> x = b; // copy
  lu.solve(x);

  cout << "x (solution) = " << x << "\n";

  // Check residual r = A*x - b
  Vector<double> r(N);
  for (size_t i=0;i<N;i++) {
    double s = 0.0;
    for (size_t j=0;j<N;j++) s += A(i,j) * x(j);
    r(i) = s - b(i);
  }
  // compute norm
  double maxabs = 0.0;
  for (size_t i=0;i<N;i++) maxabs = std::max(maxabs, std::abs(r(i)));

  cout << "max |A*x - b| = " << maxabs << "\n";

  // Also compute inverse and show A * A^{-1} approx I
  Matrix<double> Ainv = lu.inverse();
  cout << "A^{-1} =\n" << Ainv << "\n";

  Matrix<double> I = A * Ainv;
  cout << "A * A^{-1} =\n" << I << "\n";

  return 0;
}
