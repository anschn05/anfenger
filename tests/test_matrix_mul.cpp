#include <iostream>
#include <cassert>

#include "matrix.hpp"

using namespace ASC_bla;
using std::cout; using std::endl;

int main()
{
  // Simple 2x2 test
  Matrix<double> A(2,2), B(2,2);
  A(0,0)=1; A(0,1)=2; A(1,0)=3; A(1,1)=4;
  B(0,0)=5; B(0,1)=6; B(1,0)=7; B(1,1)=8;

  Matrix<double> C = A * B;

  // Expected: [1*5+2*7=19, 1*6+2*8=22; 3*5+4*7=43, 3*6+4*8=50]
  assert(C(0,0) == 19);
  assert(C(0,1) == 22);
  assert(C(1,0) == 43);
  assert(C(1,1) == 50);

  cout << "Matrix multiplication test passed." << endl;
  return 0;
}
