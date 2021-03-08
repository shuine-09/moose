#!/usr/bin/env python

import pandas
import matplotlib.pyplot as plt
import numpy as np
import math

ax = plt.gca()

data = pandas.read_csv('0.0_out.csv')
u = data['time']
F = data['resid_x_left']
plt.plot(u,-F,label='l = 20 mm without pressure',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)

data = pandas.read_csv('200_out.csv')
u = data['time']
F = data['resid_x_left']
plt.plot(u,-F,label='l = 20 mm',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)

data = pandas.read_csv('400_out.csv')
u = data['time']
F = data['resid_x_left']
plt.plot(u,-F,label='l = 10 mm',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)

data = pandas.read_csv('800_out.csv')
u = data['time']
F = data['resid_x_left']
plt.plot(u,-F,label='l = 5 mm',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)

data = pandas.read_csv('1600_out.csv')
u = data['time']
F = data['resid_x_left']
plt.plot(u,-F,label='l = 2.5 mm',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)


plt.xlabel("u[mm]")
plt.ylabel("F/A (MPa)")

plt.legend()

# plt.xlim(left=0)
# plt.ylim(bottom=0)

ax.legend(loc='best')
plt.savefig('linear_pressure_l.pdf')
