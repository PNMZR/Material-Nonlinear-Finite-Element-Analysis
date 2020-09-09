# -*- coding: utf-8 -*-
#
# Example 2.3 Two nonlinear springs (Newton-Raphson method)
#
import numpy as np

tol = 1.0e-5
iter = 0
c = 0
u = np.array([0, 0])
uold = u
f = np.array([0, 100])
P = np.array([300*u[0]**2+400*u[0]*u[1]-200*u[1]**2+150*u[0]-100*u[1],
              200*u[0]**2-400*u[0]*u[1]+200*u[1]**2-100*u[0]+100*u[1]])
R = f - P
conv = (R[0]**2+R[1]**2)/(1+f[0]**2+f[1]**2)
print('\niter u1 u2 conv c');
print('\n%3d %7.5f %7.5f %12.3e %7.5f' % (iter, u[0], u[1], conv, c))
while conv > tol and iter < 20:
    Kt = np.array([[600*u[0]+400*u[1]+150, 400*(u[0]-u[1])-100],
                   [  400*(u[0]-u[1])-100, 400*(u[1]-u[0])+100]])
    delu = np.linalg.solve(Kt, R)
    u = uold + delu
    P = np.array([300*u[0]**2+400*u[0]*u[1]-200*u[1]**2+150*u[0]-100*u[1],
                  200*u[0]**2-400*u[0]*u[1]+200*u[1]**2-100*u[0]+100*u[1]])
    R = f - P
    conv = (R[0]**2+R[1]**2)/(1+f[0]**2+f[1]**2)
    c = abs(0.9-u[1])/abs(0.9-uold[1])**2
    uold = u
    iter += 1
    print('\n%3d %7.5f %7.5f %12.3e %7.5f' % (iter, u[0], u[1], conv, c))
