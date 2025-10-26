#include <iostream>
#include "../src/matrix_inverse_neu.hpp"
using namespace ASC_bla;

int main() {
  Matrix<double> A(2, 2), B(2, 2);
  for (size_t i = 0; i < 2; i++)
    for (size_t j = 0; j < 2; j++) {
      A(i,j) = i + j;
      B(i,j) = i-j;
    }

  auto C = A * B;
  auto A_inv = Inverse(A);
  std::cout << "A^-1:\n" << A_inv << "\n\n";
  std::cout << "A:\n" << A << "\n\nB:\n" << B << "\n\nA*B:\n" << C << std::endl;
}