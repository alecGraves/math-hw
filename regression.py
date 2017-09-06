'''
regression.py:
simple linear regression algroithm using stochastic gradient descent.

author: shadySource
'''
# set 1:
# x = [10, 35, 59, 81, 106]
# y = [-3, -9, -14, -19, -25]

# set 2:
x = [10, 26, 40, 53, 69]
y = [3, -1, -4, -8, -10]

# Task specific data correction:
y = [-1*i for i in y]
y = [(i+2) for i in y]
x = [(i+2) for i in x]

a = .5

def dC():
    grad = 0
    for i in range(len(x)):
        # Derivative of slope regression function (y = ax) w/ respect to a:
        grad = grad + 2*(a*x[i] - y[i])*(a*y[i]-x[i])/(a**2 + 1)**2
    return grad

for i in range(1000):
    grad = dC()
    a = a + 1e-6*grad
    print(a, grad, i)
    input()
