# search for library like bla.cpython-312-x86_64-linux-gnu.so in the build directory:
import sys
import os

# Add the build directory to Python path
build_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'build')
sys.path.append(build_path)
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



