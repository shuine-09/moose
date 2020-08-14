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

# data = pandas.read_csv('bar_p02.csv')
# u = data['time']
# F = data['resid_x']
# plt.plot(u,F,label='pressure = 0.2',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)
#
# data = pandas.read_csv('bar_p04.csv')
# u = data['time']
# F = data['resid_x']
# plt.plot(u,F,label='pressure = 0.4',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)
#
# data = pandas.read_csv('bar_p06.csv')
# u = data['time']
# F = data['resid_x']
# plt.plot(u,F,label='pressure = 0.6',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)
#
# data = pandas.read_csv('bar_p08.csv')
# u = data['time']
# F = data['resid_x']
# plt.plot(u,F,label='pressure = 0.8',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)
#
# data = pandas.read_csv('bar_p10.csv')
# u = data['time']
# F = data['resid_x']
# plt.plot(u,F,label='pressure = 1.0',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)

# data = pandas.read_csv('bar_p08_l5.csv')
# u = data['time']
# F = data['resid_x']
# plt.plot(u,F,label='l = 5mm',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)
#
# data = pandas.read_csv('bar_p08_l10.csv')
# u = data['time']
# F = data['resid_x']
# plt.plot(u,F,label='l = 10mm',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)
#
# data = pandas.read_csv('bar_p08_l20.csv')
# u = data['time']
# F = data['resid_x']
# plt.plot(u,F,label='l = 20mm',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)
#
# data = pandas.read_csv('bar_p08_l50.csv')
# u = data['time']
# F = data['resid_x']
# plt.plot(u,F,label='l = 50mm',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)


data = pandas.read_csv('bar_p04_l50.csv')
u = data['time']
F = data['resid_x']
plt.plot(u,F,label='l = 50mm',linewidth=2.0,markerfacecolor='none',marker='>',markevery=50)

plt.xlabel("u[mm]")
plt.ylabel("F/A (MPa)")

plt.legend()

plt.xlim(left=0)
plt.ylim(bottom=0)

ax.legend(loc='best')
plt.savefig('linear_pressure_l.pdf')
