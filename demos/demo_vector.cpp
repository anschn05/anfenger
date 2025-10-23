#include <iostream>

#include "../src/vector.hpp"

namespace bla = ASC_bla;


int main()
{
  size_t n = 10;
  bla::Vector<double> x(n), y(n);

  for (size_t i = 0; i < x.size(); i++)
    {
      x(i) = i;
      y(i) = 10;
    }

  bla::Vector<double> z = x-y;
  
  std::cout << "x-y = " << z << std::endl;
}
