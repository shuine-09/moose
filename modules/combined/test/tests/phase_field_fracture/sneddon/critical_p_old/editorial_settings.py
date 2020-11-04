import matplotlib as mpl
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy.spatial import distance
from scipy import stats
from scipy.stats import t
import numpy as np

#mpl.rcParams['pdf.use14corefonts'] = True
#mpl.rcParams['ps.useafm'] = True
#mpl.rcParams['text.usetex'] = True

mpl.rcParams['axes.axisbelow'] = True
mpl.rcParams['axes.grid'] = False
mpl.rcParams['axes.titlesize'] = 13
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.edgecolor'] = '555555'
mpl.rcParams['axes.labelcolor'] = '555555'
mpl.rcParams['axes.formatter.useoffset'] = True
mpl.rcParams['axes.formatter.use_mathtext'] = True

mpl.rcParams['font.size'] = 10
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['lines.markersize'] = 8
mpl.rcParams['lines.markeredgewidth'] = 1.
#mpl.rcParams['figure.titlesize'] = 'large'
mpl.rcParams['axes.titlepad'] = 10.0     # pad between axes and title in points
#mpl.rcParams['axes.labelpad'] = 10.0
mpl.rcParams['axes.labelpad'] = 5.0

mpl.rcParams['hatch.linewidth'] = 0.5

mpl.rcParams['grid.color']= 'D4D4D4' #'b0b0b0' #whitesmoke' # 'gainsboro'
mpl.rcParams['grid.linestyle'] = ':'
mpl.rcParams['grid.linewidth'] = 0.25

mpl.rcParams['figure.figsize'] = [7.0, 4.0]
mpl.rcParams['figure.dpi'] = 80
mpl.rcParams['savefig.dpi'] = 100

#mpl.rcParams['errorbar.capsize'] = 3

mpl.rcParams['legend.frameon'] = True
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.loc'] = 'upper left'
mpl.rcParams['legend.numpoints'] = 1
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['legend.framealpha'] = None
mpl.rcParams['legend.scatterpoints'] = 1
mpl.rcParams['legend.edgecolor'] = 'inherit'
mpl.rcParams['legend.handletextpad'] = 0.05
#mpl.rcParams['legend.columnspacing'] = 0.1
mpl.rcParams['legend.columnspacing'] = 1.0
mpl.rcParams['legend.labelspacing'] = 0.2

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['xtick.color'] = '555555'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['ytick.color'] = '555555'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

mpl.rcParams['savefig.transparent'] = True

def putLegendLines(ax,gridlinecolor='silver'):
    ax.grid(True,which="major",ls="--",c=gridlinecolor,lw=0.5)
    ax.grid(True,which="minor",ls=":",c=gridlinecolor,lw=0.5)
    ax.minorticks_on()

def printStatistics(expected,predicted,message):
    print ("d= ", "%.1f"%(distance.euclidean(predicted, expected)),
           "dM= ", "%.1f"%(distance.cityblock(predicted, expected)),
           "R2= ", "%.3f"%(r2_score(predicted, expected)),
           "RMSE= ", "%.3f"%(distance.euclidean(predicted, expected)/np.sqrt(len(expected))),
           "rRMSE= ", "%.2f"%(100*distance.euclidean(1, np.divide(predicted,expected))/np.sqrt(len(expected))),
           message)


def printStatistics_v2(expected,predicted,message):
#    print ("d, dM, R2, RMSE, rRMSE, Function")
    print ("%.1f"%(distance.euclidean(predicted, expected)), " & ", # "d= ",
           "%.1f"%(distance.cityblock(predicted, expected)), " & ", # "dM= ",
           "%.3f"%(r2_score(predicted, expected)), " & ", # "R2= ",
           "%.1f"%(distance.euclidean(predicted, expected)/np.sqrt(len(expected))), " & ", # "RMSE= ",
           "%.1f"%(100*distance.euclidean(1, np.divide(predicted,expected))/np.sqrt(len(expected))), " \\\ ", # "rRMSE= ",
           message)


def make_panes_transparent(ax):
    # make the panes transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    [t.set_va('center') for t in ax.get_yticklabels()]
    [t.set_ha('left') for t in ax.get_yticklabels()]
    [t.set_va('center') for t in ax.get_xticklabels()]
    [t.set_ha('right') for t in ax.get_xticklabels()]
    [t.set_va('center') for t in ax.get_zticklabels()]
    [t.set_ha('left') for t in ax.get_zticklabels()]

    # background
#    ax.grid(False)
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    
#     tick placement
    ax.xaxis._axinfo['tick']['inward_factor'] = 0
    ax.xaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.yaxis._axinfo['tick']['inward_factor'] = 0
    ax.yaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.zaxis._axinfo['tick']['inward_factor'] = 0
    ax.zaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.zaxis._axinfo['tick']['outward_factor'] = 0.4

#    from matplotlib.ticker import MultipleLocator

#    ax.xaxis.set_major_locator(MultipleLocator(300))
#    ax.yaxis.set_major_locator(MultipleLocator(20))
#    ax.zaxis.set_major_locator(MultipleLocator(2.0))
#
    ax.view_init(elev=25, azim=-45)
