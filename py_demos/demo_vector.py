import sys, os
sys.path.append(r"C:\Users\Emanuel\Uni\Sci-Comp\anfenger\build\Debug")
try:
    from bla import Vector
except Exception as e:
    import traceback; traceback.print_exc()

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



