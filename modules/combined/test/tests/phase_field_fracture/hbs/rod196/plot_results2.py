#!/usr/bin/env python

# Wen Jiang
# Date:

import editorial_settings

import csv
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import cm
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import defaultdict
from os import system
import numpy as np
import scipy
from scipy.interpolate import griddata
from scipy.spatial import distance
from scipy import stats
from scipy.stats import t
from matplotlib import gridspec
from textwrap import wrap

from matplotlib.ticker import MultipleLocator

import pandas as pd
from shapely.geometry import LineString

cjam= plt.cm.binary(np.linspace(0.85,0.2,4))
msList= ['o','s','D','^','>','v','<','d','p','h','H','8','P','*','X','+','x']

tk_array= np.arange(300,3010,10)

benchmark = False
high_fidelity = False
correlation = False
agr567 = False
agr2 = False
fgr = False

fig = plt.figure(figsize=[6.5,5.5])
gs = gridspec.GridSpec(1,1)
ax = fig.add_subplot(gs[0])

data = pd.read_csv("bubble_pressure_r0.25_gas200_por10_ext0_rod196.csv",header=None)
time = data[0]
pressure = data[1]/1e6

x = time.to_numpy()
f = pressure.to_numpy()
x2 = x*0 + 86
g = f

first_line = LineString(np.column_stack((x, f)))
second_line = LineString(np.column_stack((x2, g)))
intersection = first_line.intersection(second_line)

if intersection.geom_type == 'MultiPoint':
    plt.plot(*LineString(intersection).xy, 'o')
elif intersection.geom_type == 'Point':
    plt.plot(*intersection.xy, 'o')

s = 'Critical pressure = ' + str(round(*intersection.xy[1], 2)) + 'MPa'

ax.text(intersection.xy[0][0]+25, intersection.xy[1][0]+0, s)

# ax.text(80, 50, 'No fracuture occurs')

ax.plot(time,pressure,'-k',c=cjam[0])

ax.set_xlabel(r"Time (second)")
ax.set_ylabel(r"Pressure (MPa)")
#ax.set_title("\n".join(wrap("Radius=0.25 $\mu$m, Initial gas pressure = 100 MPa, Porosity = 10$\%$, External pressure = 0 MPa", 50)))
ax.set_xlim(left=0)

ax.minorticks_on()
plt.savefig('bubble_pressure_r0.25_ext0_rod196.pdf', bbox_inches='tight');
plt.close(fig)

if (fgr):

    #AGR2 PIE Ag
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    labels = ['5-1-3', '5-3-3', '6-4-2', '6-4-3']
    pie = [19.3e-2, 76.7e-2, 21.8e-2, 31.2e-2]
    parfume = [12.0e-2, 23.4e-2, 12.2e-2, 9.3e-2]
    bison_ag = [1.29e-01, 2.54e-01, 1.32e-01, 1.00e-01]
    cs = [2.694086e-07,4.921005e-06,0,1.71e-11]
    sr = [1.224010e-03,2.706668e-03,1.034830e-03,7.43e-4]

    x = np.arange(len(labels))  # the label locations
    width = 0.25  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width, pie, width, label='PIE',color ='black',edgecolor = 'black')
    rects2 = ax.bar(x, parfume, width, label='Bison',color = 'silver',edgecolor = 'black')
    rects3 = ax.bar(x + width, bison_ag, width, label='PARFUME',color = 'grey',edgecolor = 'black')

    ax.set_yscale('log')
    ax.set_ylabel(r'Release Fraction')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(top = 1)
    ax.legend()

    def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate(r'${:.3f}$'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)
    autolabel(rects3)

    fig.tight_layout()

    plt.savefig('../figures/agr2_silver.pdf', bbox_inches='tight');
    plt.close(fig)

    #AGR2 PIE Cs
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    labels = ['Cesium', 'Strontium']
    pie_min = [7.26e-5, 6.35e-5]
    parfume = [4.77e-4, 1.81e-4]
    bison = [2.3816e-04,2.7e-3]
    pie_max = [1.18e-4, 3.32e-4]

    x = np.arange(len(labels))  # the label locations
    width = 0.2  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - 1.5*width, pie_min, width, label='PIE Min',color ='whitesmoke',edgecolor = 'black')
    rects2 = ax.bar(x - 0.5*width, parfume, width, label='PARFUME',color ='silver',edgecolor = 'black')
    rects3 = ax.bar(x + 0.5*width, bison, width, label='Bison',color ='grey',edgecolor = 'black')
    rects4 = ax.bar(x + 1.5*width, pie_max, width, label='PIE Max',color ='black',edgecolor = 'black')

    ax.set_yscale('log')
    ax.set_ylabel(r'Release Fraction')
    ax.set_ylim(top = 0.01)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    def as_si(x, ndp):
        s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
        m, e = s.split('e')
        return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))

    def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate(r'${}$'.format(as_si(height,2)),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)
    autolabel(rects3)
    autolabel(rects4)

    fig.tight_layout()

    plt.savefig('../figures/agr2_cs_sr.pdf', bbox_inches='tight');
    plt.close(fig)

    R = 4.2750e-04
    # AGR567 C561 Ag
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr567_561_Ag_2000.csv")
    conc = data['conc_Ag']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr567_561_Ag_0500.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,'-k',label=r"$1 \times 10^7$ s",c=cjam[3])

    data = pd.read_csv("agr567_561_Ag_1000.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,'--k',label=r"$2 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr567_561_Ag_1500.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,'-.k',label=r"$3 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr567_561_Ag_2000.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,':k',label=r"$4 \times 10^7$ s",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2485, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6140, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7778, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8655, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9532, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr567_561_silver.pdf', bbox_inches='tight');
    plt.close(fig)

    # C561 Cs
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr567_561_Cs_2000.csv")
    conc = data['conc_Cs']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr567_561_Cs_0500.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,'-k',label=r"$1 \times 10^7$ s",c=cjam[3])

    data = pd.read_csv("agr567_561_Cs_1000.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,'--k',label=r"$2 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr567_561_Cs_1500.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,'-.k',label=r"$3 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr567_561_Cs_2000.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,':k',label=r"$4 \times 10^7$ s",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2485, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6140, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7778, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8655, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9532, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr567_561_cesium.pdf', bbox_inches='tight');
    plt.close(fig)

    # C561 Sr
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr567_561_Sr_2000.csv")
    conc = data['conc_Sr']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr567_561_Sr_0500.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,'-k',label=r"$1 \times 10^7$ s",c=cjam[3])

    data = pd.read_csv("agr567_561_Sr_1000.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,'--k',label=r"$2 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr567_561_Sr_1500.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,'-.k',label=r"$3 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr567_561_Sr_2000.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,':k',label=r"$4 \times 10^7$ s",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2485, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6140, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7778, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8655, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9532, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr567_561_strontium.pdf', bbox_inches='tight');
    plt.close(fig)

    R = 4.3125e-04
    # AGR2 C513 Ag
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_513_Ag_2000.csv")
    conc = data['conc_Ag']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr2_513_Ag_0500.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,'-k',label=r"$1 \times 10^7$ s",c=cjam[3])

    data = pd.read_csv("agr2_513_Ag_1000.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,'--k',label=r"$2 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr2_513_Ag_1500.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,'-.k',label=r"$3 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr2_513_Ag_2000.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,':k',label=r"$4 \times 10^7$ s",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2474, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6094, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7709, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8586, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9497, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_513_silver.pdf', bbox_inches='tight');
    plt.close(fig)

    # C513 Cs
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_513_Cs_2000.csv")
    conc = data['conc_Cs']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr2_513_Cs_0500.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,'-k',label=r"$1 \times 10^7$ s",c=cjam[3])

    data = pd.read_csv("agr2_513_Cs_1000.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,'--k',label=r"$2 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr2_513_Cs_1500.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,'-.k',label=r"$3 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr2_513_Cs_2000.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,':k',label=r"$4 \times 10^7$ s",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2474, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6094, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7709, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8586, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9497, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_513_cesium.pdf', bbox_inches='tight');
    plt.close(fig)

    # C561 Sr
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_513_Sr_2000.csv")
    conc = data['conc_Sr']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr2_513_Sr_0500.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,'-k',label=r"$1 \times 10^7$ s",c=cjam[3])

    data = pd.read_csv("agr2_513_Sr_1000.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,'--k',label=r"$2 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr2_513_Sr_1500.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,'-.k',label=r"$3 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr2_513_Sr_2000.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,':k',label=r"$4 \times 10^7$ s",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2474, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6094, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7709, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8586, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9497, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_513_strontium.pdf', bbox_inches='tight');
    plt.close(fig)

    # AGR2 C642 Ag
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_642_Ag_2000.csv")
    conc = data['conc_Ag']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr2_642_Ag_0500.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,'-k',label=r"$1 \times 10^7$ s",c=cjam[3])

    data = pd.read_csv("agr2_642_Ag_1000.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,'--k',label=r"$2 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr2_642_Ag_1500.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,'-.k',label=r"$3 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr2_642_Ag_2000.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,':k',label=r"$4 \times 10^7$ s",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2474, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6094, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7709, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8586, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9497, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_642_silver.pdf', bbox_inches='tight');
    plt.close(fig)

    # C642 Cs
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_513_Cs_2000.csv")
    conc = data['conc_Cs']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr2_642_Cs_0500.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,'-k',label=r"$1 \times 10^7$ s",c=cjam[3])

    data = pd.read_csv("agr2_642_Cs_1000.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,'--k',label=r"$2 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr2_642_Cs_1500.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,'-.k',label=r"$3 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr2_642_Cs_2000.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,':k',label=r"$4 \times 10^7$ s",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2474, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6094, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7709, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8586, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9497, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_642_cesium.pdf', bbox_inches='tight');
    plt.close(fig)

    # C642 Sr
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_642_Cs_2000.csv")
    conc = data['conc_Cs']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr2_642_Sr_0500.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,'-k',label=r"$1 \times 10^7$ s",c=cjam[3])

    data = pd.read_csv("agr2_642_Sr_1000.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,'--k',label=r"$2 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr2_642_Sr_1500.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,'-.k',label=r"$3 \times 10^7$ s",c=cjam[2])

    data = pd.read_csv("agr2_642_Sr_2000.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,':k',label=r"$4 \times 10^7$ s",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2474, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6094, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7709, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8586, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9497, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_642_strontium.pdf', bbox_inches='tight');
    plt.close(fig)

    # AGR2 C513 failed Ag
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_513_Ag_2000.csv")
    conc = data['conc_Ag']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr2_513_fail_Ag_2000.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,'-.k',label=r"failed particle",c=cjam[1])

    data = pd.read_csv("agr2_513_Ag_2000.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,'-k',label=r"intact particle",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2474, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6094, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7709, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8586, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9497, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_513_fail_silver.pdf', bbox_inches='tight');
    plt.close(fig)

    # AGR2 C513 failed Cs
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_513_Cs_2000.csv")
    conc = data['conc_Cs']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr2_513_fail_Cs_2000.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,'-.k',label=r"failed particle",c=cjam[1])

    data = pd.read_csv("agr2_513_Cs_2000.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,'-k',label=r"intact particle",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2474, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6094, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7709, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8586, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9497, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_513_fail_cesium.pdf', bbox_inches='tight');
    plt.close(fig)

    # AGR2 C513 failed Sr
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_513_Sr_2000.csv")
    conc = data['conc_Sr']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr2_513_fail_Sr_2000.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,'-.k',label=r"failed particle",c=cjam[1])

    data = pd.read_csv("agr2_513_Sr_2000.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,'-k',label=r"intact particle",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2474, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6094, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7709, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8586, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9497, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_513_fail_strontium.pdf', bbox_inches='tight');
    plt.close(fig)

    # AGR2 C642 failed Ag
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_642_Ag_2000.csv")
    conc = data['conc_Ag']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr2_642_fail_Ag_2000.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,'-.k',label=r"failed particle",c=cjam[1])

    data = pd.read_csv("agr2_642_Ag_2000.csv")
    x = data['x']/R
    conc = data['conc_Ag']
    ax.plot(x,conc,'-k',label=r"intact particle",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2474, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6094, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7709, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8586, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9497, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_642_fail_silver.pdf', bbox_inches='tight');
    plt.close(fig)

    # AGR2 C642 failed Cs
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_642_Cs_2000.csv")
    conc = data['conc_Cs']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr2_642_fail_Cs_2000.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,'-.k',label=r"failed particle",c=cjam[1])

    data = pd.read_csv("agr2_642_Cs_2000.csv")
    x = data['x']/R
    conc = data['conc_Cs']
    ax.plot(x,conc,'-k',label=r"intact particle",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2474, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6094, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7709, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8586, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9497, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_642_fail_cesium.pdf', bbox_inches='tight');
    plt.close(fig)

    # AGR2 C642 failed Sr
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_642_Sr_2000.csv")
    conc = data['conc_Sr']
    ymax = conc.max()
    plt.plot([2.1335e-04/R, 2.1335e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.1225e-04/R, 3.1225e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.5265e-04/R, 3.5265e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([3.8785e-04/R, 3.8785e-04/R], [0, ymax], color='whitesmoke')
    plt.plot([4.3125e-04/R, 4.3125e-04/R], [0, ymax], color='whitesmoke')

    data = pd.read_csv("agr2_642_fail_Sr_2000.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,'-.k',label=r"failed particle",c=cjam[1])

    data = pd.read_csv("agr2_642_Sr_2000.csv")
    x = data['x']/R
    conc = data['conc_Sr']
    ax.plot(x,conc,'-k',label=r"intact particle",c=cjam[0])

    plt.xlabel('r / r$_0$ (-)')
    plt.ylabel('Concentration (mol/m$^3$)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=1)
    plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

    plt.text(0.2474, 1.01, 'Kernel', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.6094, 1.01, 'Buffer', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.7709, 1.01, 'IPyC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.8586, 1.01, 'SiC', horizontalalignment='center', transform=ax.transAxes)
    plt.text(0.9497, 1.01, 'OPyC', horizontalalignment='center', transform=ax.transAxes)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_642_fail_strontium.pdf', bbox_inches='tight');
    plt.close(fig)

if (agr2):
    # C513
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_513.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress_min']/1e6
    ax.plot(fluence,stress,'-k',label=r"Intact Particle",c=cjam[0])

    data = pd.read_csv("agr2_513_cracking.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'--k',label=r"IPyC Cracking",c=cjam[1])

    data = pd.read_csv("agr2_513_asphericity.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'-.k',label=r"Aspherical Particle",c=cjam[2])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_513_stress.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_513.csv")
    fluence = data['max_fluence']/1e25
    buffer_disp = data['buffer_disp']*1e3
    ipyc_disp = data['ipyc_disp']*1e3
    ax.plot(fluence,buffer_disp+0.3125,'--k',label=r"Buffer",c=cjam[1])
    ax.plot(fluence,ipyc_disp+0.3125,'-.k',label=r"IPyC",c=cjam[1])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Radius (mm)")
    ax.set_xlim(left=0)

    ax2 = ax.twinx()
    ax2.plot(fluence,(ipyc_disp-buffer_disp)*1e3,'-k',label=r'Gap',c=cjam[0])
    ax2.set_ylabel("Gap Width($\mu m$)")

    ax.plot(np.nan, '-k', label = 'Gap')
    ax.legend(loc="best")
    ax.minorticks_on()
    plt.savefig('../figures/agr2_513_gap.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_513.csv")
    fluence = data['max_fluence']/1e25
    center_line_temp = data['center_line_temp']
    outer_temp = data['outer_temp']

    ax.plot(fluence,center_line_temp-outer_temp,'-k',label=r"AGR-2 Compact 513",c=cjam[0])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Temperature Differential (K)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_513_temp.pdf', bbox_inches='tight');
    plt.close(fig)

    #533
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_533.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress_min']/1e6
    ax.plot(fluence,stress,'-k',label=r"Intact Particle",c=cjam[0])

    data = pd.read_csv("agr2_533_cracking.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'--k',label=r"IPyC Cracking",c=cjam[1])

    data = pd.read_csv("agr2_533_asphericity.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'-.k',label=r"Aspherical Particle",c=cjam[2])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_533_stress.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_533.csv")
    fluence = data['max_fluence']/1e25
    buffer_disp = data['buffer_disp']*1e3
    ipyc_disp = data['ipyc_disp']*1e3
    ax.plot(fluence,buffer_disp+0.3125,'--k',label=r"Buffer",c=cjam[1])
    ax.plot(fluence,ipyc_disp+0.3125,'-.k',label=r"IPyC",c=cjam[1])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Radius (mm)")
    ax.set_xlim(left=0)

    ax2 = ax.twinx()
    ax2.plot(fluence,(ipyc_disp-buffer_disp)*1e3,'-k',label=r'Gap',c=cjam[0])
    ax2.set_ylabel("Gap Width($\mu m$)")

    ax.plot(np.nan, '-k', label = 'Gap')
    ax.legend(loc="best")
    ax.minorticks_on()
    plt.savefig('../figures/agr2_533_gap.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_533.csv")
    fluence = data['max_fluence']/1e25
    center_line_temp = data['center_line_temp']
    outer_temp = data['outer_temp']

    ax.plot(fluence,center_line_temp-outer_temp,'-k',label=r"AGR-2 Compact 533",c=cjam[0])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Temperature Differential (K)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_533_temp.pdf', bbox_inches='tight');
    plt.close(fig)

    #642
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_642.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress_min']/1e6
    ax.plot(fluence,stress,'-k',label=r"Intact Particle",c=cjam[0])

    data = pd.read_csv("agr2_642_cracking.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'--k',label=r"IPyC Cracking",c=cjam[1])

    data = pd.read_csv("agr2_642_asphericity.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'-.k',label=r"Aspherical Particle",c=cjam[2])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_642_stress.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_642.csv")
    fluence = data['max_fluence']/1e25
    buffer_disp = data['buffer_disp']*1e3
    ipyc_disp = data['ipyc_disp']*1e3
    ax.plot(fluence,buffer_disp+0.3125,'--k',label=r"Buffer",c=cjam[1])
    ax.plot(fluence,ipyc_disp+0.3125,'-.k',label=r"IPyC",c=cjam[1])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Radius (mm)")
    ax.set_xlim(left=0)

    ax2 = ax.twinx()
    ax2.plot(fluence,(ipyc_disp-buffer_disp)*1e3,'-k',label=r'Gap',c=cjam[0])
    ax2.set_ylabel("Gap Width($\mu m$)")

    ax.plot(np.nan, '-k', label = 'Gap')
    ax.legend(loc="best")
    ax.minorticks_on()
    plt.savefig('../figures/agr2_642_gap.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_642.csv")
    fluence = data['max_fluence']/1e25
    center_line_temp = data['center_line_temp']
    outer_temp = data['outer_temp']

    ax.plot(fluence,center_line_temp-outer_temp,'-k',label=r"AGR-2 Compact 642",c=cjam[0])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Temperature Differential (K)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_642_temp.pdf', bbox_inches='tight');
    plt.close(fig)

    #643
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_643.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress_min']/1e6
    ax.plot(fluence,stress,'-k',label=r"Intact Particle",c=cjam[0])

    data = pd.read_csv("agr2_643_cracking.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'--k',label=r"IPyC Cracking",c=cjam[1])

    data = pd.read_csv("agr2_643_asphericity.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'-.k',label=r"Aspherical Particle",c=cjam[2])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_643_stress.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_643.csv")
    fluence = data['max_fluence']/1e25
    buffer_disp = data['buffer_disp']*1e3
    ipyc_disp = data['ipyc_disp']*1e3
    ax.plot(fluence,buffer_disp+0.3125,'--k',label=r"Buffer",c=cjam[1])
    ax.plot(fluence,ipyc_disp+0.3125,'-.k',label=r"IPyC",c=cjam[1])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Radius (mm)")
    ax.set_xlim(left=0)

    ax2 = ax.twinx()
    ax2.plot(fluence,(ipyc_disp-buffer_disp)*1e3,'-k',label=r'Gap',c=cjam[0])
    ax2.set_ylabel("Gap Width($\mu m$)")

    ax.plot(np.nan, '-k', label = 'Gap')
    ax.legend(loc="best")
    ax.minorticks_on()
    plt.savefig('../figures/agr2_643_gap.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr2_643.csv")
    fluence = data['max_fluence']/1e25
    center_line_temp = data['center_line_temp']
    outer_temp = data['outer_temp']

    ax.plot(fluence,center_line_temp-outer_temp,'-k',label=r"AGR-2 Compact 643",c=cjam[0])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Temperature Differential (K)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/agr2_643_temp.pdf', bbox_inches='tight');
    plt.close(fig)

if (agr567):
    # C513
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr567_513.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress_min']/1e6
    ax.plot(fluence,stress,'-k',label=r"Intact Particle",c=cjam[0])

    data = pd.read_csv("agr567_513_cracking.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'--k',label=r"IPyC Cracking",c=cjam[1])

    data = pd.read_csv("agr567_513_asphericity.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'-.k',label=r"Aspherical Particle",c=cjam[2])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/agr567_513_stress.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr567_513.csv")
    fluence = data['max_fluence']/1e25
    buffer_disp = data['buffer_disp']*1e3
    ipyc_disp = data['ipyc_disp']*1e3
    ax.plot(fluence,buffer_disp+0.3125,'--k',label=r"Buffer",c=cjam[1])
    ax.plot(fluence,ipyc_disp+0.3125,'-.k',label=r"IPyC",c=cjam[1])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Radius (mm)")
    ax.set_xlim(left=0)

    ax2 = ax.twinx()
    ax2.plot(fluence,(ipyc_disp-buffer_disp)*1e3,'-k',label=r'Gap',c=cjam[0])
    ax2.set_ylabel("Gap Width($\mu m$)")

    ax.plot(np.nan, '-k', label = 'Gap')
    ax.legend(loc="best")
    ax.minorticks_on()
    plt.savefig('../figures/agr567_513_gap.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr567_513.csv")
    fluence = data['max_fluence']/1e25
    center_line_temp = data['center_line_temp']
    outer_temp = data['outer_temp']

    ax.plot(fluence,center_line_temp-outer_temp,'-k',label=r"AGR-5/6/7 Compact 513",c=cjam[0])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Temperature Differential (K)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/agr567_513_temp.pdf', bbox_inches='tight');
    plt.close(fig)

    # C523
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr567_523.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress_min']/1e6
    ax.plot(fluence,stress,'-k',label=r"Intact Particle",c=cjam[0])

    data = pd.read_csv("agr567_523_cracking.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'--k',label=r"IPyC Cracking",c=cjam[1])

    data = pd.read_csv("agr567_523_asphericity.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'-.k',label=r"Aspherical Particle",c=cjam[2])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/agr567_523_stress.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr567_523.csv")
    fluence = data['max_fluence']/1e25
    buffer_disp = data['buffer_disp']*1e3
    ipyc_disp = data['ipyc_disp']*1e3
    ax.plot(fluence,buffer_disp+0.3125,'--k',label=r"Buffer",c=cjam[1])
    ax.plot(fluence,ipyc_disp+0.3125,'-.k',label=r"IPyC",c=cjam[1])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Radius (mm)")
    ax.set_xlim(left=0)

    ax2 = ax.twinx()
    ax2.plot(fluence,(ipyc_disp-buffer_disp)*1e3,'-k',label=r'Gap',c=cjam[0])
    ax2.set_ylabel("Gap Width($\mu m$)")

    ax.plot(np.nan, '-k', label = 'Gap')
    ax.legend(loc="best")
    ax.minorticks_on()
    plt.savefig('../figures/agr567_523_gap.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr567_523.csv")
    fluence = data['max_fluence']/1e25
    center_line_temp = data['center_line_temp']
    outer_temp = data['outer_temp']

    ax.plot(fluence,center_line_temp-outer_temp,'-k',label=r"AGR-5/6/7 Compact 523",c=cjam[0])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Temperature Differential (K)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/agr567_523_temp.pdf', bbox_inches='tight');
    plt.close(fig)

    # C561
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr567_561.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress_min']/1e6
    ax.plot(fluence,stress,'-k',label=r"Intact Particle",c=cjam[0])

    data = pd.read_csv("agr567_561_cracking.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'--k',label=r"IPyC Cracking",c=cjam[1])

    data = pd.read_csv("agr567_561_asphericity.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'-.k',label=r"Aspherical Particle",c=cjam[2])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/agr567_561_stress.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr567_561.csv")
    fluence = data['max_fluence']/1e25
    buffer_disp = data['buffer_disp']*1e3
    ipyc_disp = data['ipyc_disp']*1e3
    ax.plot(fluence,buffer_disp+0.3125,'--k',label=r"Buffer",c=cjam[1])
    ax.plot(fluence,ipyc_disp+0.3125,'-.k',label=r"IPyC",c=cjam[1])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Radius (mm)")
    ax.set_xlim(left=0)

    ax2 = ax.twinx()
    ax2.plot(fluence,(ipyc_disp-buffer_disp)*1e3,'-k',label=r'Gap',c=cjam[0])
    ax2.set_ylabel("Gap Width($\mu m$)")

    ax.plot(np.nan, '-k', label = 'Gap')
    ax.legend(loc="best")
    ax.minorticks_on()
    plt.savefig('../figures/agr567_561_gap.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("agr567_561.csv")
    fluence = data['max_fluence']/1e25
    center_line_temp = data['center_line_temp']
    outer_temp = data['outer_temp']

    ax.plot(fluence,center_line_temp-outer_temp,'-k',label=r"AGR-5/6/7 Compact 561",c=cjam[0])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Temperature Differential (K)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/agr567_561_temp.pdf', bbox_inches='tight');
    plt.close(fig)

if (correlation):
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("bison_973.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_IT']/1e6
    ax.plot(fluence,stress,'-k',label=r"Intact Particle",c=cjam[0])

    data = pd.read_csv("cracking_973.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'-.k',label=r"IPyC Cracking",c=cjam[0])

    plt.annotate("", xy=(1., 0), xytext=(1, 390),
             arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", color='silver'))
    plt.text(1, 195, r'$\bar{\sigma}_{2D}$',
         {'color': 'black', 'fontsize': 12, 'ha': 'center', 'va': 'center',
          'bbox': dict(boxstyle="round", fc="white", ec="white", pad=0.2)})

    plt.annotate("", xy=(0.9, 0), xytext=(0.9, -450),
             arrowprops=dict(arrowstyle="<->", connectionstyle="arc3",color='silver'))
    plt.text(0.9, -225, r'$\bar{\sigma}_{1D}$',
         {'color': 'black', 'fontsize': 12, 'ha': 'center', 'va': 'center',
          'bbox': dict(boxstyle="round", fc="white", ec="white", pad=0.2)})

    plt.hlines(0,0,5,colors='silver', linestyles='dotted')

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/bison_cracking_correlation.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("bison_973.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_IT']/1e6
    ax.plot(fluence,stress,'-k',label=r"Spherical Particle",c=cjam[0])

    data = pd.read_csv("asphericity_973.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'-.k',label=r"Aspherical Particle",c=cjam[0])

    plt.annotate("", xy=(1.2, 0), xytext=(1.2, -630),
             arrowprops=dict(arrowstyle="<->", connectionstyle="arc3",color='silver'))
    plt.text(1.2, -315, r'$\bar{\sigma}_{2D}$',
         {'color': 'black', 'fontsize': 12, 'ha': 'center', 'va': 'center',
          'bbox': dict(boxstyle="round", fc="white", ec="white", pad=0.2)})

    plt.annotate("", xy=(0.9, 0), xytext=(0.9, -450),
             arrowprops=dict(arrowstyle="<->", connectionstyle="arc3",color='silver'))
    plt.text(0.9, -225, r'$\bar{\sigma}_{1D}$',
         {'color': 'black', 'fontsize': 12, 'ha': 'center', 'va': 'center',
          'bbox': dict(boxstyle="round", fc="white", ec="white", pad=0.2)})

    plt.annotate("", xy=(3.4, -450), xytext=(3.4, -155),
             arrowprops=dict(arrowstyle="<->", connectionstyle="arc3",color='silver'))
    plt.text(3.4, -220, r'$\Delta\bar{\sigma}_{1D}$',
         {'color': 'black', 'fontsize': 12, 'ha': 'center', 'va': 'center',
          'bbox': dict(boxstyle="round", fc="white", ec="white", pad=0.2)})

    plt.annotate("", xy=(3.5, -630), xytext=(3.5, -280),
             arrowprops=dict(arrowstyle="<->", connectionstyle="arc3",color='silver'))
    plt.text(3.5, -550, r'$\Delta\bar{\sigma}_{2D}$',
         {'color': 'black', 'fontsize': 12, 'ha': 'center', 'va': 'center',
          'bbox': dict(boxstyle="round", fc="white", ec="white", pad=0.2)})

    plt.hlines(0,0,2,colors='silver', linestyles='dotted')
    plt.hlines(-630,1.2,4,colors='silver', linestyles='dotted')
    plt.hlines(-450,0.9,4,colors='silver', linestyles='dotted')

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/bison_asphericity_correlation.pdf', bbox_inches='tight');
    plt.close(fig)


if (high_fidelity):
    # BISON PARFUME cracking
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("bison_973.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_IT']/1e6
    ax.plot(fluence,stress,'-k',label=r"Intact 700$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1273.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_IT']/1e6
    ax.plot(fluence,stress,'--k',label=r"Intact 1000$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1573.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_IT']/1e6
    ax.plot(fluence,stress,'-.k',label=r"Intact 1300$^\circ$C",c=cjam[0])

    data = pd.read_csv("cracking_973.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'-k',label=r"IPyC Cracking 700$^\circ$C",c=cjam[2])

    data = pd.read_csv("cracking_1273.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'--k',label=r"IPyC Cracking 1000$^\circ$C",c=cjam[2])

    data = pd.read_csv("cracking_1573.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'-.k',label=r"IPyC Cracking 1300$^\circ$C",c=cjam[2])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/bison_cracking.pdf', bbox_inches='tight');
    plt.close(fig)

    # BISON PARFUME asphericity
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("bison_973.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_IT']/1e6
    ax.plot(fluence,stress,'-k',label=r"Spherical Particle 700$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1273.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_IT']/1e6
    ax.plot(fluence,stress,'--k',label=r"Spherical Particle 1000$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1573.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_IT']/1e6
    ax.plot(fluence,stress,'-.k',label=r"Spherical Particle 1300$^\circ$C",c=cjam[0])

    data = pd.read_csv("asphericity_973.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'-k',label=r"Aspherical Particle 700$^\circ$C",c=cjam[2])

    data = pd.read_csv("asphericity_1273.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'--k',label=r"Aspherical Particle 1000$^\circ$C",c=cjam[2])

    data = pd.read_csv("asphericity_1573.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_stress']/1e6
    ax.plot(fluence,stress,'-.k',label=r"Aspherical Particle 1300$^\circ$C",c=cjam[2])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/bison_asphericity.pdf', bbox_inches='tight');
    plt.close(fig)

if(benchmark):
    # BISON SIC stress
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("bison_973.csv")
    fluence = data['max_fluence']/1e25
    SiC = data['SiC_IT']/1e6
    ax.plot(fluence,SiC,'-k',label=r"700$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1273.csv")
    fluence = data['max_fluence']/1e25
    SiC = data['SiC_IT']/1e6
    ax.plot(fluence,SiC,'--k',label=r"1000$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1573.csv")
    fluence = data['max_fluence']/1e25
    SiC = data['SiC_IT']/1e6
    ax.plot(fluence,SiC,'-.k',label=r"1300$^\circ$C",c=cjam[0])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/bison_sic_stress.pdf', bbox_inches='tight');
    plt.close(fig)

    # BISON PARFUME IPYC stress
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("bison_973.csv")
    fluence = data['max_fluence']/1e25
    stress = data['IPYC_IT']/1e6
    ax.plot(fluence,stress,'-k',label=r"Bison 700$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1273.csv")
    fluence = data['max_fluence']/1e25
    stress = data['IPYC_IT']/1e6
    ax.plot(fluence,stress,'--k',label=r"Bison 1000$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1573.csv")
    fluence = data['max_fluence']/1e25
    stress = data['IPYC_IT']/1e6
    ax.plot(fluence,stress,'-.k',label=r"Bison 1300$^\circ$C",c=cjam[0])

    data = pd.read_csv('stress_973.csv',header=None, delimiter=r"\s+")
    fluence = data[1]
    stress = data[4]
    ax.plot(fluence,stress,'-k',label=r"PARFUME 700$^\circ$C",c=cjam[2])

    data = pd.read_csv('stress_1273.csv',header=None, delimiter=r"\s+")
    fluence = data[1]
    stress = data[4]
    ax.plot(fluence,stress,'--k',label=r"PARFUME 1000$^\circ$C",c=cjam[2])

    data = pd.read_csv('stress_1573.csv',header=None, delimiter=r"\s+")
    fluence = data[1]
    stress = data[4]
    ax.plot(fluence,stress,'-.k',label=r"PARFUME 1300$^\circ$C",c=cjam[2])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/bison_parfume_ipyc_stress.pdf', bbox_inches='tight');
    plt.close(fig)

    # BISON PARFUME SIC stress
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("bison_973.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_IT']/1e6
    ax.plot(fluence,stress,'-k',label=r"Bison 700$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1273.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_IT']/1e6
    ax.plot(fluence,stress,'--k',label=r"Bison 1000$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1573.csv")
    fluence = data['max_fluence']/1e25
    stress = data['SiC_IT']/1e6
    ax.plot(fluence,stress,'-.k',label=r"Bison 1300$^\circ$C",c=cjam[0])

    data = pd.read_csv('stress_973.csv',header=None, delimiter=r"\s+")
    fluence = data[1]
    stress = data[6]
    ax.plot(fluence,stress,'-k',label=r"PARFUME 700$^\circ$C",c=cjam[2])

    data = pd.read_csv('stress_1273.csv',header=None, delimiter=r"\s+")
    fluence = data[1]
    stress = data[6]
    ax.plot(fluence,stress,'--k',label=r"PARFUME 1000$^\circ$C",c=cjam[2])

    data = pd.read_csv('stress_1573.csv',header=None, delimiter=r"\s+")
    fluence = data[1]
    stress = data[6]
    ax.plot(fluence,stress,'-.k',label=r"PARFUME 1300$^\circ$C",c=cjam[2])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/bison_parfume_sic_stress.pdf', bbox_inches='tight');
    plt.close(fig)

    # BISON PARFUME OPYC stress
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("bison_973.csv")
    fluence = data['max_fluence']/1e25
    stress = data['OPYC_IT']/1e6
    ax.plot(fluence,stress,'-k',label=r"Bison 700$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1273.csv")
    fluence = data['max_fluence']/1e25
    stress = data['OPYC_IT']/1e6
    ax.plot(fluence,stress,'--k',label=r"Bison 1000$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1573.csv")
    fluence = data['max_fluence']/1e25
    stress = data['OPYC_IT']/1e6
    ax.plot(fluence,stress,'-.k',label=r"Bison 1300$^\circ$C",c=cjam[0])

    data = pd.read_csv('stress_973.csv',header=None, delimiter=r"\s+")
    fluence = data[1]
    stress = data[8]
    ax.plot(fluence,stress,'-k',label=r"PARFUME 700$^\circ$C",c=cjam[2])

    data = pd.read_csv('stress_1273.csv',header=None, delimiter=r"\s+")
    fluence = data[1]
    stress = data[8]
    ax.plot(fluence,stress,'--k',label=r"PARFUME 1000$^\circ$C",c=cjam[2])

    data = pd.read_csv('stress_1573.csv',header=None, delimiter=r"\s+")
    fluence = data[1]
    stress = data[8]
    ax.plot(fluence,stress,'-.k',label=r"PARFUME 1300$^\circ$C",c=cjam[2])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Stress (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/bison_parfume_opyc_stress.pdf', bbox_inches='tight');
    plt.close(fig)

    # BISON PARFUME pressure
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("bison_973.csv")
    fluence = data['max_fluence']/1e25
    pressure = data['gas_pressure']/1e6
    ax.plot(fluence,pressure,'-k',label=r"Bison 700$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1273.csv")
    fluence = data['max_fluence']/1e25
    pressure = data['gas_pressure']/1e6
    ax.plot(fluence,pressure,'--k',label=r"Bison 1000$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1573.csv")
    fluence = data['max_fluence']/1e25
    pressure = data['gas_pressure']/1e6
    ax.plot(fluence,pressure,'-.k',label=r"Bison 1300$^\circ$C",c=cjam[0])

    data = pd.read_csv('pressure_973.csv',header=None, delimiter=r"\s+")
    fluence = data[0]
    pressure = data[3]
    ax.plot(fluence,pressure,'-k',label=r"PARFUME 700$^\circ$C",c=cjam[2])

    data = pd.read_csv('pressure_1273.csv',header=None, delimiter=r"\s+")
    fluence = data[0]
    pressure = data[3]
    ax.plot(fluence,pressure,'--k',label=r"PARFUME 1000$^\circ$C",c=cjam[2])

    data = pd.read_csv('pressure_1573.csv',header=None, delimiter=r"\s+")
    fluence = data[0]
    pressure = data[3]
    ax.plot(fluence,pressure,'-.k',label=r"PARFUME 1300$^\circ$C",c=cjam[2])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Inner Pressure (MPa)")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/bison_parfume_pressure.pdf', bbox_inches='tight');
    plt.close(fig)

    # BISON PARFUME gas
    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    data = pd.read_csv("bison_973.csv")
    fluence = data['max_fluence']/1e25
    gas = data['fis_gas_released']
    ax.plot(fluence,gas,'-k',label=r"Bison 700$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1273.csv")
    fluence = data['max_fluence']/1e25
    gas = data['fis_gas_released']
    ax.plot(fluence,gas,'--k',label=r"Bison 1000$^\circ$C",c=cjam[0])

    data = pd.read_csv("bison_1573.csv")
    fluence = data['max_fluence']/1e25
    gas = data['fis_gas_released']
    ax.plot(fluence,gas,'-.k',label=r"Bison 1300$^\circ$C",c=cjam[0])

    data = pd.read_csv('gas_973.csv',header=None, delimiter=r"\s+")
    fluence = data[0]
    gas = data[4]
    ax.plot(fluence,gas,'-k',label=r"PARFUME 700$^\circ$C",c=cjam[2])

    data = pd.read_csv('gas_1273.csv',header=None, delimiter=r"\s+")
    fluence = data[0]
    gas = data[4]
    ax.plot(fluence,gas,'--k',label=r"PARFUME 1000$^\circ$C",c=cjam[2])

    data = pd.read_csv('gas_1573.csv',header=None, delimiter=r"\s+")
    fluence = data[0]
    pressure = data[4]
    ax.plot(fluence,pressure,'-.k',label=r"PARFUME 1300$^\circ$C",c=cjam[2])

    ax.set_xlabel(r"Fluence (10$^{25}$n/m$^2$)")
    ax.set_ylabel(r"Gas Moles")
    ax.legend(loc="best")
    ax.set_xlim(left=0)

    ax.minorticks_on()
    plt.savefig('../figures/bison_parfume_gas.pdf', bbox_inches='tight');
    plt.close(fig)
