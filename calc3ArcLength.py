'''
Simple math utility script
ex. usage:
    from mymath import *
    v = 8*t*N.i + (t**(3/2)+1)*N.j
    arc(v, 1, 8)
'''
from __future__ import division
from sympy import *
x, y, z, t = symbols('x y z t')

from sympy.vector import CoordSysCartesian
N = CoordSysCartesian('N')

def arc(v,i,f):
    abs_r_prime = diff(v, t).magnitude()
    print(abs_r_prime)
    int_abs_r_prime = integrate(abs_r_prime, t)
    print(int_abs_r_prime)
    def_int_abs_r_prime=integrate(abs_r_prime, (t, i, f))
    print(def_int_abs_r_prime)

v = 8*t*N.i + (t**6+2)/(t**2)*N.j
arc(v, 1, 8)
