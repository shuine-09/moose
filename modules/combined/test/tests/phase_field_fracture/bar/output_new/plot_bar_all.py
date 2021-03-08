#!/usr/bin/env python

import pandas
import matplotlib.pyplot as plt
import numpy as np
import math

ax = plt.gca()

data = pandas.read_csv('0.0_out.csv')
u = data['time']
F = data['resid_x_left']
plt.plot(u,-F,label='pressure = 0.0 MPa',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)

data = pandas.read_csv('0.2_out.csv')
u = data['time']
F = data['resid_x_left']
plt.plot(u,-F,label='pressure = 0.2 MPa',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)

data = pandas.read_csv('0.4_out.csv')
u = data['time']
F = data['resid_x_left']
plt.plot(u,-F,label='pressure = 0.4 MPa',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)

data = pandas.read_csv('0.6_out.csv')
u = data['time']
F = data['resid_x_left']
plt.plot(u,-F,label='pressure = 0.6 MPa',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)

data = pandas.read_csv('0.8_out.csv')
u = data['time']
F = data['resid_x_left']
plt.plot(u,-F,label='pressure = 0.8 MPa',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)

data = pandas.read_csv('1.0_out.csv')
u = data['time']
F = data['resid_x_left']
plt.plot(u,-F,label='pressure = 1.0 MPa',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)

plt.xlabel("u[mm]")
plt.ylabel("F/A (MPa)")

plt.legend()

# plt.xlim(left=0)
# plt.ylim(bottom=0)

ax.legend(loc='best')
plt.savefig('linear_pressure_l.pdf')
