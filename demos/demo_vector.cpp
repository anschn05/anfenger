#include <iostream>

#include "../src/vector.hpp"

namespace bla = ASC_bla;


int main()
{
  size_t n = 5;
  bla::Vector<double> x(n), y(n);

  for (size_t i = 0; i < x.Size(); i++)
    {
      x(i) = i;
      y(i) = 10;
    }

  bla::Vector<double> z1 = x-y;
  bla::Vector<double> z2 = x+y;
  
  std::cout << "x-y = " << z1 << std::endl;
  std::cout << "x+y = " << z2 << std::endl;
  
  std::cout << "type of (x+3*y) is  " << typeid(x+3*y).name() << std::endl;


  std::cout << "sizeof(x+3*y) = " << sizeof(x+3*y) << std::endl;
  
  std::cout << "<x,x> = " << dot(x,x) << std::endl;

  x.range(2,9) = 3;
  x.slice(1,5) = 10;
  
  std::cout << "x = " << x << std::endl;  
}
