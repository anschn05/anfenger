# Using pybind to apply C++ libraries also in Python

There is also the possibility to use these libraries in python. 
This happens with **pybind**.

The file, which manages the "translation" is *bind_bla.cpp* in 'src'

So if you want to use the library in a python code, you need to install a few things and follow the instructions of [Python bindings](https://jschoeberl.github.io/IntroSC/basiclinalg/bla-python2.html)

and then run the following commands:

    cmake -S . -B build -G "Visual Studio 17 2022"
    cd build
    cmake --build .


To implement the library in your code, use something like this - *according to your path of course*:

    # search for libraray like bla.cpython-312-darwin.so in the build directory:
    import sys
    sys.path.append(r"C:\Users\ansch\OneDrive\Dokumente\Uni\SciComp\anfenger\build\Debug")
    from bla import Vector

    # import from the installed ASCsoft package:
    #from ASCsoft.bla import Vector


An example code, would then look like this:

    # search for libraray like bla.cpython-312-darwin.so in the build directory:
    import sys
    sys.path.append(r"C:\Users\ansch\OneDrive\Dokumente\Uni\SciComp\anfenger\build\Debug")
    from bla import Vector

    # import from the installed ASCsoft package:
    #from ASCsoft.bla import Vector

    x = Vector(3)
    y = Vector(3)

    for i in range(len(x)):
        x[i] = i
    y[:] = 2    

    print ("x =", x)
    print ("y =", y)
    print ("x+3*y =", x+3*y)


    x = Vector(10)
    x[0:] = 1
    print (x)

    x[3:7] = 2
    print (x)

    x[0:10:2] = 3
    print (x)



