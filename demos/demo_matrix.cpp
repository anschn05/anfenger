#include <iostream>
#include "../src/matrix.hpp"
using namespace ASC_bla;

int main() {
  Matrix<double> A(2, 2), B(2, 2);
  for (size_t i = 0; i < 2; i++)
    for (size_t j = 0; j < 3; j++) {
      A(i,j) = i + j;
      B(i,j) = i-j;
    }

  auto C = A * B;
  std::cout << "A:\n" << A << "\n\nB:\n" << B << "\n\nA*B:\n" << C << std::endl;
}