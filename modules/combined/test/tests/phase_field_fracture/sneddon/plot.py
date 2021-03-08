import pandas
import matplotlib.pyplot as plt
import numpy as np
import math
import glob
import os

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

from matplotlib.ticker import MultipleLocator

cjam= plt.cm.binary(np.linspace(0.85,0.2,4))
msList= ['o','s','D','^','>','v','<','d','p','h','H','8','P','*','X','+','x']

tk_array= np.arange(300,3010,10)

dir_path = os.path.dirname(os.path.realpath(__file__))

all_files = glob.glob(dir_path + "/critical_l004/*.csv")

p = [0] * len(all_files)
length = [0] * len(all_files)
j = 0
for filename in all_files:
    df = pandas.read_csv(filename, index_col=None, header=0)
    time = df['time']
    c = df['integrated_c']
    x = filename.split("/")
    y = x[12].split("_")
    length[j] = float(y[0])
    vec = df.to_numpy()
    diff_old = 100
    for i in range(1, vec.shape[0]):
        diff_new = np.abs(vec[i,1] - vec[i-1,1])
        if diff_new > diff_old * 5 :
            p[j] = vec[i-1,0]
            j = j+1
            break
        else:
            diff_old = diff_new

a = np.linspace(0.025,0.8,100) # 100 linearly spaced numbers
th_p = np.sqrt(1.0*0.001/(1-0.2*0.2)/np.pi/a)
sigma0 = np.sqrt(0.001/0.04)

fig = plt.figure(figsize=[6.5,5.5])
gs = gridspec.GridSpec(1,1)
ax = fig.add_subplot(gs[0])

length_a = np.array(length)
p_a = np.array(p)

ax.plot(a/0.04,th_p/sigma0,'-k',c=cjam[0], label='Analytical solution')
ax.plot(length_a/0.04,p_a/sigma0,linestyle = 'None',color='tab:blue', marker='^', markersize = 10, label = 'Phase field fracture simulation l = 0.04')



all_files = glob.glob(dir_path + "/critical_l002/*.csv")

p = [0] * len(all_files)
length = [0] * len(all_files)
j = 0
for filename in all_files:
    df = pandas.read_csv(filename, index_col=None, header=0)
    time = df['time']
    c = df['integrated_c']
    x = filename.split("/")
    y = x[12].split("_")
    length[j] = float(y[0])
    vec = df.to_numpy()
    diff_old = 100
    for i in range(1, vec.shape[0]):
        diff_new = np.abs(vec[i,1] - vec[i-1,1])
        if diff_new > diff_old * 15 :
            p[j] = vec[i-1,0]
            j = j+1
            break
        else:
            diff_old = diff_new

a = np.linspace(0.025,0.8,100) # 100 linearly spaced numbers
th_p = np.sqrt(1.0*0.001/(1-0.2*0.2)/np.pi/a)
sigma0 = np.sqrt(0.001/0.02)

length_a = np.array(length)
p_a = np.array(p)

#ax.plot(a/0.02,th_p/sigma0,'-k',c=cjam[0], label='Analytical solution')
ax.plot(length_a/0.02,p_a/sigma0,linestyle = 'None',color='tab:red', marker='^', markersize = 10, label = 'Phase field fracture simulation l = 0.02')

all_files = glob.glob(dir_path + "/critical_l002_elem/*.csv")

p = [0] * len(all_files)
length = [0] * len(all_files)
j = 0
for filename in all_files:
    df = pandas.read_csv(filename, index_col=None, header=0)
    time = df['time']
    c = df['integrated_c']
    x = filename.split("/")
    y = x[12].split("_")
    length[j] = float(y[0])
    vec = df.to_numpy()
    diff_old = 100
    for i in range(1, vec.shape[0]):
        diff_new = np.abs(vec[i,1] - vec[i-1,1])
        if diff_new > diff_old * 15 :
            p[j] = vec[i-1,0]
            j = j+1
            break
        else:
            diff_old = diff_new

a = np.linspace(0.025,0.8,100) # 100 linearly spaced numbers
th_p = np.sqrt(1.0*0.001/(1-0.2*0.2)/np.pi/a)
sigma0 = np.sqrt(0.001/0.02)

length_a = np.array(length)
p_a = np.array(p)

#ax.plot(a/0.02,th_p/sigma0,'-k',c=cjam[0], label='Analytical solution')
ax.plot(length_a/0.02,p_a/sigma0,linestyle = 'None',color='tab:orange', marker='^', markersize = 10, label = 'Phase field fracture simulation l = 0.02 elem ic')

all_files = glob.glob(dir_path + "/critical_l001_elem/*.csv")

p = [0] * len(all_files)
length = [0] * len(all_files)
j = 0
for filename in all_files:
    df = pandas.read_csv(filename, index_col=None, header=0)
    time = df['time']
    c = df['integrated_c']
    x = filename.split("/")
    y = x[12].split("_")
    length[j] = float(y[0])
    vec = df.to_numpy()
    diff_old = 100
    for i in range(1, vec.shape[0]):
        diff_new = np.abs(vec[i,1] - vec[i-1,1])
        if diff_new > diff_old * 15 :
            p[j] = vec[i-1,0]
            j = j+1
            break
        else:
            diff_old = diff_new

a = np.linspace(0.025,0.8,100) # 100 linearly spaced numbers
th_p = np.sqrt(1.0*0.001/(1-0.2*0.2)/np.pi/a)
sigma0 = np.sqrt(0.001/0.01)

length_a = np.array(length)
p_a = np.array(p)

#ax.plot(a/0.02,th_p/sigma0,'-k',c=cjam[0], label='Analytical solution')
ax.plot(length_a/0.01,p_a/sigma0,linestyle = 'None',color='tab:brown', marker='^', markersize = 10, label = 'Phase field fracture simulation l = 0.01 elem ic')

all_files = glob.glob(dir_path + "/critical_l004_elem/*.csv")

p = [0] * len(all_files)
length = [0] * len(all_files)
j = 0
for filename in all_files:
    df = pandas.read_csv(filename, index_col=None, header=0)
    time = df['time']
    c = df['integrated_c']
    x = filename.split("/")
    y = x[12].split("_")
    length[j] = float(y[0])
    vec = df.to_numpy()
    diff_old = 100
    for i in range(1, vec.shape[0]):
        diff_new = np.abs(vec[i,1] - vec[i-1,1])
        if diff_new > diff_old * 15 :
            p[j] = vec[i-1,0]
            j = j+1
            break
        else:
            diff_old = diff_new

a = np.linspace(0.025,0.8,100) # 100 linearly spaced numbers
th_p = np.sqrt(1.0*0.001/(1-0.2*0.2)/np.pi/a)
sigma0 = np.sqrt(0.001/0.04)

length_a = np.array(length)
p_a = np.array(p)

#ax.plot(a/0.02,th_p/sigma0,'-k',c=cjam[0], label='Analytical solution')
ax.plot(length_a/0.04,p_a/sigma0,linestyle = 'None',color='tab:purple', marker='^', markersize = 10, label = 'Phase field fracture simulation l = 0.04 elem ic')


all_files = glob.glob(dir_path + "/critical_l001/*.csv")

p = [0] * len(all_files)
length = [0] * len(all_files)
j = 0
for filename in all_files:
    df = pandas.read_csv(filename, index_col=None, header=0)
    time = df['time']
    c = df['integrated_c']
    x = filename.split("/")
    y = x[12].split("_")
    length[j] = float(y[0])
    vec = df.to_numpy()
    diff_old = 100
    for i in range(1, vec.shape[0]):
        diff_new = np.abs(vec[i,1] - vec[i-1,1])
        if diff_new > diff_old * 15 :
            p[j] = vec[i-1,0]
            j = j+1
            break
        else:
            diff_old = diff_new

a = np.linspace(0.025,0.8,100) # 100 linearly spaced numbers
th_p = np.sqrt(1.0*0.001/(1-0.2*0.2)/np.pi/a)
sigma0 = np.sqrt(0.001/0.01)

length_a = np.array(length)
p_a = np.array(p)

#ax.plot(a/0.02,th_p/sigma0,'-k',c=cjam[0], label='Analytical solution')
ax.plot(length_a/0.01,p_a/sigma0,linestyle = 'None',color='tab:green', marker='^', markersize = 10, label = 'Phase field fracture simulation l = 0.01')



plt.xlabel('Normalized Crack Length a/$l_0$')
plt.ylabel('Normalized Pressure p/$\sigma_0$')
ax.set_ylim(bottom=0)
ax.set_xlim(left=0,right=20)
plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

plt.savefig('lefm_pf.pdf')
plt.close()

# for filename in all_files:
#     df = pandas.read_csv(filename, index_col=None, header=0)
#     time = df['time']
#     c = df['integrated_c']
#     ax = plt.gca()
#     ax.plot(time,c,'-k',c=cjam[0])
#     plt.savefig(filename+'.pdf')
#     plt.close()
