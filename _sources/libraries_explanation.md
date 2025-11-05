# How to use the libraries and functions

### Right now there are the following header-files:
- matrix.hpp
- matrixexpr.hpp
- vecexpr.hpp
- vector.hpp

The most important functions are:

for Vectors:

    x = Vector(3)           # creates Vector with length 3
    y = Vector(3)
    y[:] = 2                # y = [2,2,2]
    x[1] = 0
    x[2:3] = 1              # x = [0,1,1] 

    x+y                     # adds two vectors
    x-y                     # substracts two vectors
    print("x =",x)          # prints: x = 1,2,3
    y[:] = 2                # changes y to [2,2,2]

    
