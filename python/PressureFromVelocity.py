#!/usr/bin/env python3

rho = 1.293 # kg/m3
cd = 0.8
pi = 3.14159

A = 0.2 * 0.15 # m2
a = pi * (3/1000.0)**2 # m2
v = 1.0 # m/s

P = 0.5 * rho * (A**2) * (v**2) / (cd * a)**2 # N/m2
F = P * A

print('{} [kPa]'.format(P/1000))
print('{} [PSI]'.format(P/6895))
