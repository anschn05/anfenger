# search for libraray like bla.cpython-312-darwin.so in the build directory:
# import sys
# sys.path.append('/Users/joachim/texjs/lva/ws2324/ScientificComputing/ASC-bla/build')
# from bla import Vector

# import from the installed ASCsoft package:
#from ASCsoft.bla import Vector

import sys,os
sys.path.append(r"C:\Users\ansch\OneDrive - TU Wien\Documents\Uni\5_Semester\SciComp\anfenger\ASC-bla\build\Debug")
os.add_dll_directory(r"C:/Program Files (x86)/Microsoft Visual Studio/2022/BuildTools/VC/Tools/MSVC/14.44.35207/bin")
from bla import Vector


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



