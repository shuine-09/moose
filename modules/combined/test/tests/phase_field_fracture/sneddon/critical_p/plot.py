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
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy.spatial import distance
from scipy import stats
from scipy.stats import t
from matplotlib import gridspec

from matplotlib.ticker import MultipleLocator

cjam= plt.cm.binary(np.linspace(0.85,0.2,4))
msList= ['o','s','D','^','>','v','<','d','p','h','H','8','P','*','X','+','x']

tk_array= np.arange(300,3010,10)

dir_path = os.path.dirname(os.path.realpath(__file__))

all_files = glob.glob(dir_path + "/*.csv")

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

# p = np.sort(p)
# p = p[::-1]

#length = [0.1,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.3,0.4,0.5,0.6]
# p = [0.054,0.049,0.046,0.044,0.041,0.04,0.038,0.036,0.035,0.034,0.033,0.029,0.026,0.024]

a = np.linspace(0.1,0.8,100) # 100 linearly spaced numbers
th_p = np.sqrt(1.0*0.001/(1-0.2*0.2)/np.pi/a)

fig = plt.figure(figsize=[6.5,5.5])
gs = gridspec.GridSpec(1,1)
ax = fig.add_subplot(gs[0])

ax.plot(a,th_p,'-k',c=cjam[0], label='Analytical solution')
ax.plot(length,p,linestyle = 'None',color='tab:blue', marker='^', markersize = 10, label = 'Phase field fracture simulation')

plt.xlabel('Crack Length')
plt.ylabel('Critical Pressure')
#ax.set_ylim(bottom=0)
#ax.set_xlim(left=0,right=0.8)
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
