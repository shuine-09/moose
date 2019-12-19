#!/usr/bin/env python

import pandas
import matplotlib.pyplot as plt
import numpy as np
import math

ax = plt.gca()

CP = pandas.read_csv('moving_oxide_CP_weak_mechanics_out.csv')
time_cp = CP['time']
thickness_cp = CP['cut_data_x']

C4 = pandas.read_csv('moving_oxide_C4_weak_out.csv')
time_c4 = C4['time']
thickness_c4 = C4['cut_data_x']

plt.plot(time_cp,thickness_cp,label='CP (BISON XFEM)')
plt.plot(time_c4,thickness_c4,label='C4 (BISON XFEM)')

#plt.title("Plastic strain")
plt.xlabel("time[s]")
plt.ylabel("oxide thickness [$\mu$ m]")

ax.legend()
#ax.set_xlim(left=0)
#ax.set_ylim(bottom=0,top = 0.5)
plt.savefig('C4_and_CP_model.pdf')
