#!/usr/bin/env python

import pandas
import matplotlib.pyplot as plt
import numpy as np
import math

ax = plt.gca()

case1 = pandas.read_csv('notch_c_w0.csv')
strain1 = case1['time']
stress1 = case1['stress_yy']

case2 = pandas.read_csv('notch_c_w10.csv')
strain2 = case2['time']
stress2 = case2['stress_yy']

case3 = pandas.read_csv('notch_c_w20.csv')
strain3 = case3['time']
stress3 = case3['stress_yy']

case4 = pandas.read_csv('notch_c_w80.csv')
strain4 = case4['time']
stress4 = case4['stress_yy']

plt.plot(strain1,stress1,label='W0 = 0')
plt.plot(strain2,stress2, label='W0 = 1e7')
plt.plot(strain3,stress3, label='W0 = 2e7')
plt.plot(strain4,stress4, label='W0 = 8e7')

plt.title("Stress")
plt.xlabel("Total Strain")
plt.ylabel("Stress")

ax.legend()
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
plt.savefig('Fig11.pdf')
