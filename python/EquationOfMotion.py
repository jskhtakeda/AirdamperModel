#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import math
# np.seterr(divide='ignore', invalid='ignore')

####################

PI = 3.14159
rho_0 = 1.293 # kg/m3
p_0 = 101325 # Pa = N/m2
g = 9.8 # m/s2
gamma = 1.4

c_o = 0.8
A_o = 0.2 * 0.15 # m2

A_p = PI * (5/1000.0)**2 # m2
h_p = 0.2 # m

m = 8.0 # kg

C = (c_o * A_o) / (A_p * h_p) * math.sqrt(2 * p_0 / rho_0)
D = (m * g) / (p_0 * A_p)
print("C: {}, D:{}".format(C, D))

N = 100
t_0 = 0.0 # sec
t_1 = 0.3 # sec
dt = (t_1-t_0) / N

x = 0 # m
v = 1.0 # m/s
a = 9.8 # m/s2

Q = A_p * h_p / 0.1
Re = Q * (6 / 1000) / (1.8 * 1E-5) / A_p
print("Re: {}".format(Re)) # Re < 10 で層流

####################

def dxdt(t, x, v, a):
    print("dxdt: {}".format(v))
    return v

def dvdt(t, x, v, a):
    print("dvdt: {}".format(a))
    return a

def dadt(t, x, v, a):
    bunshi = C * np.sign(1 - a / g) * math.sqrt(D * abs(1 - a / g)) - (D * (1 - a / g) + 1) ** (1 / gamma) * v / h_p
    bunbo = (1 / gamma) * (1 - x / h_p) * (D * (1 - a / g) + 1) ** (1 / gamma - 1) * D / g
    print("dadt: {}".format(bunshi / bunbo))
    return bunshi / bunbo

t_array = np.arange(t_0, t_1, dt)
x_array = []
v_array = []
a_array = []

for t in t_array:
    x_array.append(x)
    v_array.append(v)
    a_array.append(a)
    
    k1x = dxdt(t, x, v, a)
    k1v = dvdt(t, x, v, a)
    k1a = dadt(t, x, v, a)

    k2x = dxdt(t+dt/2, x+dt/2*k1x, v+dt/2*k1v, a+dt/2*k1a)
    k2v = dvdt(t+dt/2, x+dt/2*k1x, v+dt/2*k1v, a+dt/2*k1a)
    k2a = dadt(t+dt/2, x+dt/2*k1x, v+dt/2*k1v, a+dt/2*k1a)
    
    k3x = dxdt(t+dt/2, x+dt/2*k2x, v+dt/2*k2v, a+dt/2*k2a)
    k3v = dvdt(t+dt/2, x+dt/2*k2x, v+dt/2*k2v, a+dt/2*k2a)
    k3a = dadt(t+dt/2, x+dt/2*k2x, v+dt/2*k2v, a+dt/2*k2a)

    k4x = dxdt(t+dt, x+dt*k3x, v+dt*k3v, a+dt*k3a)
    k4v = dvdt(t+dt, x+dt*k3x, v+dt*k3v, a+dt*k3a)
    k4a = dadt(t+dt, x+dt*k3x, v+dt*k3v, a+dt*k3a)
    
    x += dt / 6 * (k1x + 2 * k2x + 2* k3x + k4x)
    v += dt / 6 * (k1v + 2 * k2v + 2* k3v + k4v)
    a += dt / 6 * (k1a + 2 * k2a + 2* k3a + k4a)

    # x += dt * (k1x)
    # v += dt * (k1v)
    # a += dt * (k1a)
    
plt.plot (t_array, x_array)
plt.xlabel("t")
plt.ylabel("x")

# plt.show()
