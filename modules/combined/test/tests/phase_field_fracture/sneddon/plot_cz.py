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

all_files = glob.glob(dir_path + "/critical_p_0p05_cz/*.csv")

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
        if diff_new > diff_old * 10 :
            p[j] = vec[i-1,0]
            j = j+1
            break
        else:
            diff_old = diff_new

a = np.linspace(0.025,0.8,100) # 100 linearly spaced numbers
th_p = np.sqrt(1.0*0.001*1.12433979929/(1-0.2*0.2)/np.pi/a)
sigma0 = np.sqrt(0.001*1.12433979929/0.02)

fig = plt.figure(figsize=[6.5,5.5])
gs = gridspec.GridSpec(1,1)
ax = fig.add_subplot(gs[0])

length_a = np.array(length)
p_a = np.array(p)

#p_a=[0.0557,0.053,0.050,0.047,0.045,0.043,0.0415,0.0396,0.038,0.036668,0.0355] #old
p_a=[0.077,0.067,0.061,0.0546964,0.0512002,0.0481857,0.0461729,0.0435323,0.0416261,0.0398452,0.0381489]
length_a = [0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3]

ax.plot(a/0.02,th_p/sigma0,'-k',c=cjam[0], label='Analytical solution')

#th_p = np.sqrt(1.0*0.001/(1-0.2*0.2)/np.pi/a)
#sigma0 = np.sqrt(0.001/0.02)

#ax.plot(a/0.02,th_p/sigma0,'--k',c=cjam[0], label='Analytical solution')

ax.plot(np.divide(length_a,0.02),np.divide(p_a,sigma0),linestyle = 'None',color='tab:blue', marker='^', markersize = 10, label = 'Fracture strength $\sigma_0 = 0.05$')


#p_a=[0.08013,0.0686,0.0621889,0.0562415,0.0524696,0.0490998,0.0464172,0.0446839,0.04195,0.0404,0.039]
p_a=[0.086,0.074,0.066,0.0584047,0.0545853,0.050928,0.0483716,0.0456052,0.0432649,0.0412515,0.0394953]
length_a = [0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3]

ax.plot(np.divide(length_a,0.02),np.divide(p_a,sigma0),linestyle = 'None',color='tab:red', marker='^', markersize = 10, label = 'Fracture strength $\sigma_0 = 0.10$')


# p_a=[0.08697,0.0742,0.0663,0.0587,0.0547,0.0511,0.0478,0.0460,0.04315,0.041668,0.3989]
p_a=[0.092,0.077,0.068,0.0611858,0.0554316,0.0517115,0.0487381,0.0469872,0.0439204,0.0418765,0.0403929]
length_a = [0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3]

ax.plot(np.divide(length_a,0.02),np.divide(p_a,sigma0),linestyle = 'None',color='tab:green', marker='^', markersize = 10, label = 'Fracture strength $\sigma_0 = 0.15$')

# p_a=[0.092,0.077,0.068,0.06025,0.05585,0.05249,0.0491046,0.0462,0.0439,0.0418765,0.0400937]
p_a=[0.095,0.079,0.069,0.0611858,0.0562779,0.052495,0.0492878,0.0469872,0.0442482,0.0418765,0.0406921]
length_a = [0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3]

ax.plot(np.divide(length_a,0.02),np.divide(p_a,sigma0),linestyle = 'None',color='tab:orange', marker='^', markersize = 10, label = 'Fracture strength $\sigma_0 = 0.20$')


# all_files = glob.glob(dir_path + "/critical_p_0p1_cz/*.csv")
#
# p = [0] * len(all_files)
# length = [0] * len(all_files)
# j = 0
# for filename in all_files:
#     df = pandas.read_csv(filename, index_col=None, header=0)
#     time = df['time']
#     c = df['integrated_c']
#     x = filename.split("/")
#     y = x[12].split("_")
#     length[j] = float(y[0])
#     vec = df.to_numpy()
#     diff_old = 100
#     for i in range(1, vec.shape[0]):
#         diff_new = np.abs(vec[i,1] - vec[i-1,1])
#         if diff_new > diff_old * 2.5 :
#             p[j] = vec[i-1,0]
#             j = j+1
#             break
#         else:
#             diff_old = diff_new
#
# a = np.linspace(0.025,0.8,100) # 100 linearly spaced numbers
# th_p = np.sqrt(1.0*0.001/(1-0.2*0.2)/np.pi/a)
# sigma0 = np.sqrt(0.001/0.02)
#
# length_a = np.array(length)
# p_a = np.array(p)
#
# #ax.plot(a/0.02,th_p/sigma0,'-k',c=cjam[0], label='Analytical solution')
# ax.plot(length_a/0.02,p_a/sigma0,linestyle = 'None',color='tab:red', marker='^', markersize = 10, label = 'Fracture strength $\sigma_0 = 0.1$')
#
# all_files = glob.glob(dir_path + "/critical_p_0p2_cz/*.csv")
#
# p = [0] * len(all_files)
# length = [0] * len(all_files)
# j = 0
# for filename in all_files:
#     df = pandas.read_csv(filename, index_col=None, header=0)
#     time = df['time']
#     c = df['integrated_c']
#     x = filename.split("/")
#     y = x[12].split("_")
#     length[j] = float(y[0])
#     vec = df.to_numpy()
#     diff_old = 100
#     for i in range(1, vec.shape[0]):
#         diff_new = np.abs(vec[i,1] - vec[i-1,1])
#         if diff_new > diff_old * 2.5 :
#             p[j] = vec[i-1,0]
#             j = j+1
#             break
#         else:
#             diff_old = diff_new
#
# a = np.linspace(0.025,0.8,100) # 100 linearly spaced numbers
# th_p = np.sqrt(1.0*0.001/(1-0.2*0.2)/np.pi/a)
# sigma0 = np.sqrt(0.001/0.02)
#
# length_a = np.array(length)
# p_a = np.array(p)
#
# #ax.plot(a/0.02,th_p/sigma0,'-k',c=cjam[0], label='Analytical solution')
# ax.plot(length_a/0.02,p_a/sigma0,linestyle = 'None',color='tab:green', marker='^', markersize = 10, label = 'Fracture strength $\sigma_0 = 0.2$')


plt.xlabel('Normalized Crack Length a/$l_0$')
plt.ylabel('Normalized Pressure p/$\sigma_0$')
ax.set_ylim(bottom=0,top=0.75)
ax.set_xlim(left=0,right=20)
plt.legend(framealpha=1.0, loc="upper right", bbox_to_anchor=(1, 1))

plt.savefig('lefm_pf_cz.pdf')
plt.close()

# for filename in all_files:
#     df = pandas.read_csv(filename, index_col=None, header=0)
#     time = df['time']
#     c = df['integrated_c']
#     ax = plt.gca()
#     ax.plot(time,c,'-k',c=cjam[0])
#     plt.savefig(filename+'.pdf')
#     plt.close()
