'''
Simple math utility script
ex. usage:
    from calc3 import *
    v = 8*t*N.i + (t**(3/2)+1)*N.j
    arc(v, 1, 8)
'''
from __future__ import division
from sympy import *
x, y, z, t = symbols('x y z t')
u, v, w = symbols('u v w')

from sympy.vector import CoordSysCartesian
N = CoordSysCartesian('N')

def arc(v,i,f):
    abs_r_prime = diff(v, t).magnitude()
    print(abs_r_prime)
    int_abs_r_prime = integrate(abs_r_prime, t)
    print(int_abs_r_prime)
    def_int_abs_r_prime=integrate(abs_r_prime, (t, i, f))
    print(def_int_abs_r_prime)

def grad3d(f, point=None, display=False):
    '''returns the 3D gradient'''
    grad = diff(f,x)*N.i + diff(f,y)*N.j + diff(f,z)*N.k
    if point is not None:
        x_0, y_0, z_0 = point
        grad = grad.subs([(x, x_0), (y, y_0), (z, z_0)])
    if display:
        print(grad)
    return grad

def tangent_plane(f,point):
    '''
    f in implicit form (f(x,y,z)=0)
    returns tangent plane formula at specific point
    (output=0 is the equation for the tangent plane)
    '''
    x_0, y_0, z_0 = point
    grad_f = grad3d(f, point=point)
    return grad_f.dot((x-x_0)*N.i+(y-y_0)*N.j+(z-z_0)*N.k)

def directional_diff(f, direction, point=None, display=False):
    '''returns directional derivitive of f'''
    direction = direction[0]*N.i + direction[1]*N.j + direction[2]*N.k
    direction = direction.normalize()
    diff = direction.dot(grad3d(f))
    if point is not None:
        x_0, y_0, z_0 = point
        diff = diff.subs([(x, x_0), (y, y_0), (z, z_0)])
    if display:
        print(diff)
    return diff

def partials(f, display=False):
    '''returns partial derivitives of f'''
    partials = {}
    for symbol in f.free_symbols:
        deriv = diff(f, symbol)
        if display:
            print("df/d"+str(symbol)+' =', deriv)
        partials.update({str(symbol): deriv})
    return partials

def critical_points(f):
    '''prints critical plints'''
    # where the gradient = 0:
    print(solve([
            Eq(diff(f,x), 0),
            Eq(diff(f,y), 0),
            Eq(diff(f,z), 0)]))
    # discontinuities
    print(solve([
            Eq(diff(f,x), zoo),
            Eq(diff(f,y), zoo),
            Eq(diff(f,z), zoo)]))
    # points = solve(Eq(grad3d(f), 0), x)
    # where gradient DNE:
    # points.append(Eq(grad3d(f), zoo))

def double_integral(f, i, f1, method=1):
    '''
    prints the single and double integral of f from i to f1

    ex. double_integral(x + y, x**2, sqrt(125*x))
    '''
    if method is 1:
        first_var = list(f.free_symbols)
        second_var = list(i.free_symbols)[0]
        try:
            first_var.remove(list(i.free_symbols)[0])
        except ValueError:
            pass
        intersecs = solve(Eq(i, f1))
        print(intersecs)
        i1 = integrate(f, (first_var, i, f1))
        i2 = integrate(i1, (second_var, intersecs[0], intersecs[1]))
        print(i1)
        print(i2)

def surface_area(f, dom1=None, dom2=None, mode=0):
    '''
    prints the surface area of 3d surface
    ex. 
    '''
    if mode is 0:
        print('integrating:', sqrt(diff(f,x)**2 + diff(f,y)**2 + 1))
        print(integrate(sqrt(diff(f,x)**2 + diff(f,y)**2 + 1),
                        (x, dom1[0], dom1[1]),
                        (y, dom2[0], dom2[1])))
    elif mode is 1:
        r, a = symbols('r a')
        radius = sqrt(solve(f, x**2 + y**2)[0])
        f = diff(f,x)**2 + diff(f,y)**2+1
        coeff = -1*solve(f, (y**2))[0].subs(x, 0)**(-1)
        igrl = integrate(sqrt(coeff*r**2+1)*r, (r, 0, radius), (a, 0, 2*pi))
        print(igrl)
    elif mode is 2:
        r = Symbol('r')
        f = (sqrt(diff(f,x)**2 + diff(f,y)**2 + 1)*r).subs(x**2 + y**2, r**2)
        print(f)
        print(integrate(f,(r, dom1[0], dom1[1]),(t, dom2[0], dom2[1])))


def grad_and_direc_dir(f, p, q):
    '''
    ind the gradient 
    ∇f(x, y)
    of f at the given point P. 
    (b) Use the gradient to find the directional derivative 
    Duf(x, y)
    of f at P in the direction from P to Q. 
    '''
    p = list(p)
    direction = [q[i]-p for i, p in enumerate(p)]
    while len(direction) < 3:
        direction.append(0)
        p.append(0)
    grad3d(f, p, True)
    print(directional_diff(f, direction).subs([(x, p[0]),
                                (y, p[1]), (z, p[2])]))

def max_dir_and_value(f, point):
    '''
     Find the direction for which the directional derivative 
        of the function f is a max
     Find the max value of the directional derivative
    '''
    grad = grad3d(f, point)
    # max dir:
    grad_dir = grad.normalize()
    print(grad_dir)
    #Max value:
    # directional_diff(f,
    #     [grad_dir.dot(N.i), grad_dir.dot(N.j), grad_dir.dot(N.k)],
    #     point,
    #     display=True)
    print(grad.magnitude())


# v = 8*t*N.i + (t**6+2)/(t**2)*N.j
# arc(v, 1, 8)

# f = x**2 + 1**2 - z
# pt = (2,1,5)
# print(tangent_plane(f, pt))

# print(diff(6*exp(x)*ln(y), x).subs([(x,0), (y,E)]))

# f = ln(sqrt(x**2 + y**2))
# print(directional_diff(f, (1,0,0)).subs([(x,8),(y,15)]))
# print(directional_diff(f, (0,1,0)).subs([(x,8),(y,15)]))

# x = u*t
# y = exp(u + 5*v + 3*w +7*t)
# z = u + 1/2*v + 4*t
# f = x + 2*y**2 - z**2
# print(diff(f,u))
# print(diff(f,v))
# print(diff(f,w))
# print(diff(f,t))

# # dz/dx and dz/dy
# f = exp(4*y*z)*ln(x) + y*exp(9*x*z)-y*z
# print(-1*diff(f, x)/diff(f, z))
# print(-1*diff(f, y)/diff(f, z))

# w = (7*x)**(9*y + 5*z)
# partials(w)
# f = sin(8*x)*cos(4*y + z*2)
# grad_and_direc_dir(f, (1,1,1), (2,-1,0))

# f = sqrt(81*x**2 + 81*y**2)
# max_dir_and_value(f, (3,4,0))

# f = 2*x**2 + y**2 -z
# p = (-2,4,24)
# print(tangent_plane(f, p))

# reimann sum example
# f = 8*x*(4-y)
# summ = 0
# for i in range(0,3):
#     for j in range(1,4):
#         summ += f.subs([(x,i),(y,j)])
# print(summ)


# print(integrate(1, (y, 0, 1/sqrt(x-1)), (x, 2, 5)))


# double_integral(y**2,2-x, x**2)


# print(integrate(sqrt(25-x**2), (x, 0, 3), (y, -1, 4)))

########################################################
# Find the surface area described.
# the part of the surface 
# z = xy
#  in the first octant that lies within the cylinder 
# x2 + y2 = a2,
# a > 0:
#######################################################
# a, r = symbols('a r')
# z = x*y
# surface_area(z, (0,a), (0,pi/2), 2)

# f = x*y
# f = integrate(f, (z, 0, x*y), (y, 0, -x+1), (x, 0, 1))

# print(f)
