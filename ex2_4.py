# -*- coding: utf-8 -*-
#
# Example 2.4 Divergence of the Newton-Raphson method.
#
import numpy as np
import matplotlib.pyplot as plt


xdata = np.zeros(41)
ydata = np.zeros(41)
tol = 1.0e-5
iter = 0
u = 0.5
uold = u
c = 0
P = u + np.arctan(5*u)
R = -P
conv= R**2
xdata[0] = u
ydata[0] = P
while conv > tol and iter < 20:  # 最大迭代次数为20次
    Kt = 1+5*(np.cos(np.arctan(5*u)))**2
    delu = R/Kt
    u = uold + delu
    P = u + np.arctan(5*u)
    R = -P
    conv = R**2
    uold = u
    iter = iter + 1;
    xdata[2*iter-1] = u
    ydata[2*iter-1] = 0
    xdata[2*iter] = u
    ydata[2*iter] = P

fig, ax = plt.subplots()
ax.plot(xdata, ydata, label='Newton–Raphson method')
x = np.arange(-1, 1.1, 0.1)
y = x + np.arctan(5*x)
y1 = np.zeros_like(x)
ax.plot(x, y, label=r"$P(u) = u + tan^{-1}(5u)$")
ax.plot(x, y1, 'b--')
ax.set_xlim(-1, 1)
ax.legend()
ax.set_xlabel("u")
ax.set_ylabel("P(u)")
ax.set_title("Divergence of the Newton-Raphson method")