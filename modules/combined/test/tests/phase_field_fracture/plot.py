#!/usr/bin/env python

import pandas
import matplotlib.pyplot as plt
import numpy as np
import math

ax = plt.gca()

# data = pandas.read_csv('b5.csv')
# u = data['time']
# F = data['resid_x']
# plt.plot(u,F,label='linear softening (l=5mm)',linewidth=2.0,markerfacecolor='none',marker='o',markevery=50)
#
# data = pandas.read_csv('b10.csv')
# u = data['time']
# F = data['resid_x']
# plt.plot(u,F,label='linear softening (l=10mm)',linewidth=2.0,markerfacecolor='none',marker='D',markevery=50)
#
# data = pandas.read_csv('b20.csv')
# u = data['time']
# F = data['resid_x']
# plt.plot(u,F,label='linear softening (l=20mm)',linewidth=2.0,markerfacecolor='none',marker='^',markevery=50)
#
# data = pandas.read_csv('b50.csv')
# u = data['time']
# F = data['resid_x']
# plt.plot(u,F,label='linear softening (l=50mm)',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)
#
# plt.xlabel("u[mm]")
# plt.ylabel("F/A (MPa)")
#
# plt.legend()
#
# plt.xlim(left=0)
# plt.ylim(bottom=0)
#
# ax.legend(loc='best')
# plt.savefig('linear.pdf')

data = pandas.read_csv('b5s.csv')
u = data['time']
F = data['resid_x']
plt.plot(u,F,label='standard phase field fracture (l=5mm)',linewidth=2.0,markerfacecolor='none',marker='o',markevery=50)

data = pandas.read_csv('b10s.csv')
u = data['time']
F = data['resid_x']
plt.plot(u,F,label='standard phase field fracture (l=10mm)',linewidth=2.0,markerfacecolor='none',marker='D',markevery=50)

data = pandas.read_csv('b20s.csv')
u = data['time']
F = data['resid_x']
plt.plot(u,F,label='standard phase field fracture (l=20mm)',linewidth=2.0,markerfacecolor='none',marker='^',markevery=50)

data = pandas.read_csv('b50s.csv')
u = data['time']
F = data['resid_x']
plt.plot(u,F,label='lstandard phase field fracture (l=50mm)',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)

plt.xlabel("u[mm]")
plt.ylabel("F/A (MPa)")

plt.legend()

plt.xlim(left=0)
plt.ylim(bottom=0)

ax.legend(loc='best')
plt.savefig('standard.pdf')
