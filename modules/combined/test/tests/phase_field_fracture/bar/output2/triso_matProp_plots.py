#!/usr/bin/env python

# Aysenur Toptan
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
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy.spatial import distance
from scipy import stats
from scipy.stats import t
from matplotlib import gridspec

from matplotlib.ticker import MultipleLocator

import pandas as pd

cjam= plt.cm.binary(np.linspace(0.85,0.2,4))
msList= ['o','s','D','^','>','v','<','d','p','h','H','8','P','*','X','+','x']

tk_array= np.arange(300,3010,10)

bufferThermal= False
ucoThermal= False
sicThermal= False
graphiteThermal= False


kernelMech= False
bufferMech= True
pycMech= False


if(bufferThermal):

    def matProp_k(rho,kinit=0.5,rhoinit=1000,ktheo=4.0,rhotheo=2250):
        return kinit*ktheo*rhotheo*(rhotheo-rhoinit) / ( ktheo*rhotheo*(rhotheo-rho) + kinit * rho * (rho-rhoinit))

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    ax.set_xlim([900,1000])
    #    ax.set_ylim([0,10])

    rho_array= np.arange(900,1010,10)
    ax.plot(rho_array,matProp_k(rho_array),'-k',label=r"thermal conductivity")

#    ax.xaxis.set_major_locator(MultipleLocator(300))

    #editorial_settings.putLegendLines(ax)

    ax.minorticks_on()

    ax.set_xlabel(r"$\rho$ (kg/m$^3$)")
    ax.set_ylabel(r"$k$ (W/m-K)")

    ax.legend(loc="upper left")

    plt.savefig('../figures/plot_buffer_thermal.pdf', bbox_inches='tight');
    plt.close(fig)


if(ucoThermal):

    def matProp_k(tk):
        k= []
        for i in range(len(tk)):
            tc= tk[i]-273.15
            if tc<1650: k.append( 0.0132*np.exp(0.00188*tc)+4040/(464+tc) )
            else: k.append( 0.0132*np.exp(0.00188*tc)+1.9 )
        return k

    def matProp_cp(tk):
        return 52.1743 + 87.951*tk/1000 - 84.2411 * np.power((tk/1000),2) + 31.542 * np.power((tk/1000),3) - 2.6334 * np.power((tk/1000),4) - 0.71391 / np.power((tk/1000),2)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    axR = ax.twinx()

    ax.set_xlim([300,3000])
    #    ax.set_ylim([0,10])

    ax.plot(tk_array,matProp_k(tk_array),'-k',label=r"thermal conductivity")
    axR.plot(tk_array,matProp_cp(tk_array),'--k',label=r"specific heat capacity")

    ax.xaxis.set_major_locator(MultipleLocator(300))

    #editorial_settings.putLegendLines(ax)

    ax.minorticks_on()

    ax.set_xlabel(r"$T$ (K)")
    ax.set_ylabel(r"$k$ (W/m-K)")
    axR.set_ylabel(r"$C_P$ (J/mol-K)")

    ax.legend(loc="upper left")
    axR.legend(loc="upper right")

    plt.savefig('../figures/plot_uco_thermal.pdf', bbox_inches='tight');
    plt.close(fig)



if(sicThermal):

    def matProp_k(tk): return 17885/tk + 2.0

    def matProp_cp(tk): return 925.65 + 0.3772 * tk - 7.9259e-5*np.power(tk,2) - 3.1946e7 /np.power(tk,2)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    axR = ax.twinx()

    ax.set_xlim([300,3000])
    #    ax.set_ylim([0,10])

    ax.plot(tk_array,matProp_k(tk_array),'-k',label=r"thermal conductivity")
    axR.plot(tk_array,matProp_cp(tk_array),'--k',label=r"specific heat capacity")

    ax.xaxis.set_major_locator(MultipleLocator(300))

    #editorial_settings.putLegendLines(ax)

    ax.minorticks_on()

    ax.set_xlabel(r"$T$ (K)")
    ax.set_ylabel(r"$k$ (W/m-K)")
    axR.set_ylabel(r"$c_P$ (J/kg-K)")

    ax.legend(loc="upper left")
    axR.legend(loc="center right")

    plt.savefig('../figures/plot_sic_thermal.pdf', bbox_inches='tight');
    plt.close(fig)


if(graphiteThermal):
    def matProp_kunirr(tk,k100,alpha,delta): return k100 * (1-alpha*(tk-100)*np.exp(delta*tk))
    def matProp_kphi(tk,phi): return 1.0 - (0.94-0.604*(tk-273.15)/1000) * (1-np.exp(-(2.96-1.955*(tk-273.15)/1000))*phi/1.52) - (0.043*(tk-273.15)/1000-0.008*np.power((tk-273.15)/1000,8)) * phi/1.52
    def matProp_kpf(pf,tk,k100,alpha,delta):
        kp= 4.13
        beta= (kp - matProp_kunirr(tk,k100,alpha,delta)) / (kp+2*matProp_kunirr(tk,k100,alpha,delta))
        return (1+2*beta*pf+ (2*np.power(beta,3)-0.1*beta)*np.power(pf,2) + 0.05*np.power(pf,3)*np.exp(4.5*beta)) / (1-beta*pf)

    def matProp_k(tk,k100,alpha,delta):
        k= []
        for i in range(len(tk)):
            kunirr= matProp_kunirr(tk[i],k100,alpha,delta)
            kphi= 1.0
            kpf= 1.0
            krho= 1.0
            k.append(kunirr * kphi * kpf * krho)

        return k

    def matProp_cp(tk,rho=1.7e3):
        cp= []
        for i in range(len(tk)):
            tc= tk[i]-273.15
            cp.append( 1.75e6*(0.645+3.14*tc/1000 -2.809*np.power((tc/1000),2)+0.959*np.power((tc/1000),3))/rho )
        return cp

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

#    axR = ax.twinx()

    ax.set_xlim([300,3000])
    #    ax.set_ylim([0,10])

#    ax.plot(tk_array,matProp_k(tk_array,50.8,1.1819e-3,-7.8453e-4),'-k',c=cjam[0],label=r"Eq.(3.6) A3-3")
#    ax.plot(tk_array,matProp_k(tk_array,64.6,1.4079e-3,-9.0739e-4),'-k',c=cjam[2])#,label=r"Eq.(3.6) A3-3 at 1950$^\circ$C")
#    ax.plot(tk_array,matProp_k(tk_array,47.4,9.7556e-4,-6.0360e-4),'-.k',c=cjam[0],label=r"Eq.(3.6) A3-27")
#    ax.plot(tk_array,matProp_k(tk_array,64.6,1.4621e-3,-9.6050e-4),'-.k',c=cjam[2])#,label=r"Eq.(3.6) A3-27 at 1950$^\circ$C")
    ax.plot(tk_array,matProp_cp(tk_array),'-k',label=r"Eq.(3.13)")

    ax.xaxis.set_major_locator(MultipleLocator(300))

    #editorial_settings.putLegendLines(ax)

    ax.minorticks_on()

    ax.set_xlabel(r"$T$ (K)")
#    ax.set_ylabel(r"$k$ (W/m-K)")
    ax.set_ylabel(r"$c_P$ (J/kg-K)")

#    ax.legend(loc="upper left")
    ax.legend(loc="upper right")

    plt.savefig('../figures/plot_graphite_thermal.pdf', bbox_inches='tight');
    plt.close(fig)



    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])


    ax.set_xlim([300,3000])

    ax.plot(tk_array,matProp_kunirr(tk_array,50.8,1.1819e-3,-7.8453e-4),'-k',c=cjam[0],label=r"Eq.(3.7) A3-3")
    ax.plot(tk_array,matProp_kunirr(tk_array,64.6,1.4079e-3,-9.0739e-4),'-k',c=cjam[2])
    ax.plot(tk_array,matProp_kunirr(tk_array,47.4,9.7556e-4,-6.0360e-4),'-k',dashes=[8,2],c=cjam[0],label=r"Eq.(3.7) A3-27")
    ax.plot(tk_array,matProp_kunirr(tk_array,64.6,1.4621e-3,-9.6050e-4),'-k',dashes=[8,2],c=cjam[2])

    ax.xaxis.set_major_locator(MultipleLocator(300))

    ax.text(950,33,r"1950$^\circ$C")
    ax.text(1050,26,r"1800$^\circ$C")

    ax.minorticks_on()

    ax.set_xlabel(r"$T$ (K)")
    ax.set_ylabel(r"$k_{unirr}$ (W/m-K)")

    ax.legend(loc="upper left")

    plt.savefig('../figures/plot_graphite_kunirr_thermal.pdf', bbox_inches='tight');
    plt.close(fig)


    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])


    ax.set_xlim([0,4])
    phi_array= np.arange(0,4.01,0.01)

#    kphi= []
#    for i in range(len(phi_array)):
#        kphi.append( matProp_kphi(1000,phi_array[i]) )

    ax.plot(phi_array,matProp_kphi(2000,phi_array),'-k',c=cjam[0],label=r"Eq.(3.8)")
    ax.plot(phi_array,matProp_kphi(1500,phi_array),'-k',c=cjam[1])#,label=r"Eq.(3.8)")
    ax.plot(phi_array,matProp_kphi(1000,phi_array),'-k',c=cjam[2])#,label=r"Eq.(3.8)")
    ax.plot(phi_array,matProp_kphi(500,phi_array),'-k',c=cjam[3])#,label=r"Eq.(3.8)")

#    ax.xaxis.set_major_locator(MultipleLocator(300))
    ax.minorticks_on()

    ax.text(0.5,1.28,r"2000K",rotation=24)
    ax.text(0.5,0.85,r"1500K",rotation=8)
    ax.text(0.5,0.55,r"1000K",rotation=6)
    ax.text(0.5,0.24,r"500K",rotation=5)

    ax.set_xlabel(r"$\phi$ (10$^{25}$ n/m$^2$)")
    ax.set_ylabel(r"$\kappa_{\phi}$ (-)")

    ax.legend(loc="upper left")

    plt.savefig('../figures/plot_graphite_kphi_thermal.pdf', bbox_inches='tight');
    plt.close(fig)


    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])


    ax.set_xlim([0,1])
    pf_array= np.arange(0,1.01,0.01)

#    ax.plot([0,1],[0,0],'-k',c='dimgray',lw=0.75)

    ax.plot(pf_array,matProp_kpf(pf_array,2000,50.8,1.1819e-3,-7.8453e-4),'-k',c=cjam[0],label=r"Eq.(3.9) A3-3 1800$^\circ$C")
    ax.plot(pf_array,matProp_kpf(pf_array,2000,64.6,1.4079e-3,-9.0739e-4),'-k',c=cjam[1],label=r"Eq.(3.9) A3-3 1950$^\circ$C")
    ax.plot(pf_array,matProp_kpf(pf_array,2000,47.4,9.7556e-4,-6.0360e-4),'-k',c=cjam[2],label=r"Eq.(3.9) A3-27 1800$^\circ$C")
    ax.plot(pf_array,matProp_kpf(pf_array,2000,64.6,1.4621e-3,-9.6050e-4),'-k',c=cjam[3],label=r"Eq.(3.9) A3-27 1950$^\circ$C")

#    ax.xaxis.set_major_locator(MultipleLocator(300))
    ax.minorticks_on()


    ax.set_xlabel(r"$PF$ (-)")
    ax.set_ylabel(r"$\kappa_{PF}$ (-)")

    ax.legend(loc="lower left")

    plt.savefig('../figures/plot_graphite_kpf_thermal.pdf', bbox_inches='tight');
    plt.close(fig)


    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])


    ax.set_xlim([1000,2400])
    rho_array= np.arange(1000,2410,10)

#    ax.plot([0,1],[0,0],'-k',c='dimgray',lw=0.75)

    ax.plot(rho_array,rho_array/1700,'-k',c=cjam[0],label=r"Eq.(3.12)")

#    ax.xaxis.set_major_locator(MultipleLocator(300))
    ax.minorticks_on()


    ax.set_xlabel(r"$\rho$ (kg/m$^3$)")
    ax.set_ylabel(r"$\kappa_{\rho}$ (-)")

    ax.legend(loc="upper left")

    plt.savefig('../figures/plot_graphite_krho_thermal.pdf', bbox_inches='tight');
    plt.close(fig)


if(kernelMech):

    cjam= plt.cm.binary(np.linspace(0.85,0.2,6))

    def matProp_E(tk,rho,rhoth):
        tc= tk-273.15
        return 219 * (1-1e-4*tc -2.1e-7*np.power(tc,2)+3.1e-10*np.power(tc,3)-1.6e-13*np.power(tc,4))*(1.92*rho-0.92*rhoth)/rhoth

    def matProp_mu(rho,rhoth):
        return 1.35 * (1.92*rho-0.92*rhoth) / (1.66*rho-0.66*rhoth) -1

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    axR = fig.add_subplot(gs[0])

    axR.set_xlim([8.00,13.00])
    axR.set_ylim([0.2,0.4])
    rho_array= np.arange(0.85*11.0,11.1,0.10)
    axR.plot(rho_array,matProp_mu(rho_array,11.0),'-k',c=cjam[3],label=r"$\rho_{th}=11$g/cm$^3$")
    rho_array= np.arange(0.85*13.0,13.1,0.10)
    axR.plot(rho_array,matProp_mu(rho_array,13.0),'-k',c=cjam[5],label=r"$\rho_{th}=13$g/cm$^3$")

#    ax.xaxis.set_major_locator(MultipleLocator(300))
    axR.minorticks_on()

    axR.set_xlabel(r"$\rho$ (g/cm$^3$)")
    axR.set_ylabel(r"$\nu$ (-)")

    axR.legend(loc="upper left")

    plt.savefig('../figures/plot_kernel_mech_v2.pdf', bbox_inches='tight');
    plt.close(fig)


    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    ax.set_xlim([8.00,13.00])
    ax.set_ylim([0,250])
    rho_array= np.arange(0.85*11.0,11.1,0.10)
    ax.plot(rho_array,matProp_E(1000,rho_array,11.0),'-k',c=cjam[0],label=r"$\rho_{th}=11$g/cm$^3$, $T=1000K$")
    ax.plot(rho_array,matProp_E(2000,rho_array,11.0),'--k', dashes=[7,2],c=cjam[0],label=r"$\rho_{th}=11$g/cm$^3$, $T=2000K$")
    rho_array= np.arange(0.85*13.0,13.1,0.10)
    ax.plot(rho_array,matProp_E(1000,rho_array,13.0),'-k',c=cjam[3],label=r"$\rho_{th}=13$g/cm$^3$, $T=1000K$")
    ax.plot(rho_array,matProp_E(2000,rho_array,13.0),'--k', dashes=[7,2],c=cjam[3],label=r"$\rho_{th}=13$g/cm$^3$, $T=2000K$")

    ax.minorticks_on()

    ax.set_xlabel(r"$\rho$ (g/cm$^3$)")
    ax.set_ylabel(r"$E$ (GPa)")

    ax.legend(loc="upper left")

    plt.savefig('../figures/plot_kernel_mech_v1.pdf', bbox_inches='tight');
    plt.close(fig)


if(bufferMech):
    def matProp_cte(tk): return 5*(1+0.11*((tk-273.15)-400)/700)
    def matProp_eiso(phi,a1,a2,a3,a4): return a1*phi + a2*np.power(phi,2) + a3*np.power(phi,3) + a4*np.power(phi,4)
    def matProp_ks(tk,rho): return 2.0*(1+2.38*(1.9-rho))*(2.19310e-29-4.8510e-32*(tk-273.15)+4.014710e-35*np.power((tk-273.15),2))

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    ax.set_xlim([300,3000])
    ax.xaxis.set_major_locator(MultipleLocator(300))

    ax.plot(tk_array,matProp_cte(tk_array),'-k',c=cjam[0],label=r"thermal expansion coefficient")

#    ax.xaxis.set_major_locator(MultipleLocator(300))
    ax.minorticks_on()

    ax.set_xlabel(r"$T$ (K)")
    ax.set_ylabel(r"$\alpha$ (10$^{-6}$/K)")

    ax.legend(loc="upper left")

    plt.savefig('../figures/plot_buffer_mech_v2.pdf', bbox_inches='tight');
    plt.close(fig)


    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    phi_array= np.arange(0,4.01,0.01)

    ax.set_xlim([0,4])

    ax.plot([0,4],[0,0],'-k',lw=0.75,c='dimgray')
    ax.plot(phi_array,matProp_eiso(phi_array,-1.42840,-0.19563,0.18991,-0.02591),'-k',c=cjam[0],label=r"1350$^\circ$C")
    ax.plot(phi_array,matProp_eiso(phi_array,-1.52390,0.13048,0.06299,-0.01072),'--k',c=cjam[1],label=r"1032$^\circ$C")
    ax.plot(phi_array,matProp_eiso(phi_array,-1.24080,0.00175,0.08533,-0.01253),'-.k',c=cjam[2],label=r"600$^\circ$C")

#    ax.xaxis.set_major_locator(MultipleLocator(300))
    ax.minorticks_on()

    ax.set_xlabel(r"$\phi$ (10$^{25}$ n/m$^2$)")
    ax.set_ylabel(r"$\epsilon_{iso}$ (%)")

    ax.legend(loc="lower left")

    plt.savefig('../figures/plot_buffer_mech_v3.pdf', bbox_inches='tight');
    plt.close(fig)


    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    ax.set_xlim([300,3000])
    ax.xaxis.set_major_locator(MultipleLocator(300))

#    ax.plot([0,4],[0,0],'-k',lw=0.75,c='dimgray')
    ax.plot(tk_array,matProp_ks(tk_array,1.7),'-k',c=cjam[2],label=r"$\rho=1.7$g/cm$^3$")
    ax.plot(tk_array,matProp_ks(tk_array,1.8),'-k',c=cjam[1],label=r"$\rho=1.8$g/cm$^3$")
    ax.plot(tk_array,matProp_ks(tk_array,1.9),'-k',c=cjam[0],label=r"$\rho=1.9$g/cm$^3$")
    ax.plot(tk_array,matProp_ks(tk_array,2.0),'-k',c=cjam[1],label=r"$\rho=2.0$g/cm$^3$")
    ax.plot(tk_array,matProp_ks(tk_array,2.1),'-k',c=cjam[2],label=r"$\rho=2.1$g/cm$^3$")

#    ax.xaxis.set_major_locator(MultipleLocator(300))
    ax.minorticks_on()

    ax.text(2500,3.7e-28,r"$1.7$g/cm$^3$",rotation=60)
    ax.text(2500,3.1e-28,r"$1.8$g/cm$^3$",rotation=55)
    ax.text(2500,2.5e-28,r"$1.9$g/cm$^3$",rotation=50)
    ax.text(2500,1.9e-28,r"$2.0$g/cm$^3$",rotation=40)
    ax.text(2500,1.3e-28,r"$2.1$g/cm$^3$",rotation=30)

    ax.set_xlabel(r"$T$ (K)")
    ax.set_ylabel(r"$K_S$ (m$^2$/MPa-n)")

    #ax.legend(loc="upper left")

    plt.savefig('../figures/plot_buffer_mech_v4.pdf', bbox_inches='tight');
    plt.close(fig)


    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    ax.set_xlim([300,3000])
    ax.xaxis.set_major_locator(MultipleLocator(300))

#    ax.plot([0,4],[0,0],'-k',lw=0.75,c='dimgray')
#    ax.plot(tk_array,matProp_ks(tk_array,1.7),'-k',c=cjam[2])#,label=r"Eq.(3.21) $\rho=1.8$g/cm$^3$")

#    ax.xaxis.set_major_locator(MultipleLocator(300))
    ax.minorticks_on()


    ax.set_xlabel(r"$T$ (K)")
    ax.set_ylabel(r"$E$ (GPa)")

    ax.legend(loc="upper left")

    plt.savefig('../figures/plot_buffer_mech_v1.pdf', bbox_inches='tight');
    plt.close(fig)


if(pycMech):

    phi_array= np.arange(0,4.01,0.01)
    baf_array= [1.0000, 1.0212, 1.0488, 1.0769, 1.1746, 1.2787]
    baf_tan_array= [1.0000, 1.0303, 1.0769, 1.1250, 1.2258, 1.3333]

    def matProp_eiso(phi,a1,a2,a3,a4): return a1*phi + a2*np.power(phi,2) + a3*np.power(phi,3) + a4*np.power(phi,4)

    def matProp_eradial_t600(phi,baf):
        if(baf==baf_array[0]): return matProp_eiso(phi,-1.24080,0.00175,0.08533,-0.01253)
        elif(baf==baf_array[1]): return matProp_eiso(phi,-1.10640,-0.03128,0.09184,-0.01220)
        elif(baf==baf_array[2]): return matProp_eiso(phi,-0.94333,-0.03589,0.08184,-0.00958)
        elif(baf==baf_array[3]): return matProp_eiso(phi,-0.78045,-0.02975,0.06655,-0.00626)
        elif(baf==baf_array[4]): return matProp_eiso(phi,-0.15714,-0.14889,0.07546,-0.00293)
        elif(baf==baf_array[5]): return matProp_eiso(phi,0.40265,-0.16501,0.03676,0.00706)

    def matProp_eradial_t1032(phi,baf):
        if(baf==baf_array[0]): return matProp_eiso(phi,-1.52390,0.13048,0.06299,-0.01072)
        elif(baf==baf_array[1]): return matProp_eiso(phi,-2.07520,1.37845,-0.48993,0.06602)
        elif(baf==baf_array[2]): return matProp_eiso(phi,-2.00470,1.30380,-0.37280,0.04538)
        elif(baf==baf_array[3]): return matProp_eiso(phi,-1.81690,1.10850,-0.23868,0.02484)
        elif(baf==baf_array[4]): return matProp_eiso(phi,-1.18540,0.64995,0.01380,-0.01284)
        elif(baf==baf_array[5]): return matProp_eiso(phi,-0.45900,0.51172,-0.03245,-0.00142)

    def matProp_eradial_t1350(phi,baf):
        if(baf==baf_array[0]): return matProp_eiso(phi,-1.42840,-0.19563,0.18991,-0.02591)
        elif(baf==baf_array[1]): return matProp_eiso(phi,-1.54330,0.59804,-0.09997,0.00978)
        elif(baf==baf_array[2]): return matProp_eiso(phi,-1.49640,1.16621,-0.30106,0.03475)
        elif(baf==baf_array[3]): return matProp_eiso(phi,-0.89522,0.80331,-0.09009,0.00467)
        elif(baf==baf_array[4]): return matProp_eiso(phi,1.20930,-0.53861,0.43114,-0.05590)
        elif(baf==baf_array[5]): return matProp_eiso(phi,3.71620,-2.70420,1.17990,-0.13910)


    def matProp_etangential_t600(phi,baf):
        if(baf==baf_tan_array[0]): return matProp_eiso(phi,-1.24080,0.00175,0.08533,-0.01253)
        elif(baf==baf_tan_array[1]): return matProp_eiso(phi,-1.38550,0.05307,0.07620,-0.01245)
        elif(baf==baf_tan_array[2]): return matProp_eiso(phi,-1.46790,-0.02836,0.12139,-0.01948)
        elif(baf==baf_tan_array[3]): return matProp_eiso(phi,-1.64660,0.03928,0.10067,-0.01764)
        elif(baf==baf_tan_array[4]): return matProp_eiso(phi,-1.84990,-0.09358,0.18119,-0.03036)
        elif(baf==baf_tan_array[5]): return matProp_eiso(phi,-2.19190,0.02675,0.15352,-0.02972)

    def matProp_etangential_t1032(phi,baf):
        if(baf==baf_tan_array[0]): return matProp_eiso(phi,-1.52390,0.13048,0.06299,-0.01072)
        elif(baf==baf_tan_array[1]): return matProp_eiso(phi,-1.57590,0.09019,0.05306,-0.00815)
        elif(baf==baf_tan_array[2]): return matProp_eiso(phi,-1.32200,-0.51928,0.27603,-0.03465)
        elif(baf==baf_tan_array[3]): return matProp_eiso(phi,-1.18700,-0.90635,0.41046,-0.05067)
        elif(baf==baf_tan_array[4]): return matProp_eiso(phi,-0.96963,-1.59110,0.64689,-0.07682)
        elif(baf==baf_tan_array[5]): return matProp_eiso(phi,-0.81239,-2.20760,0.88496,-0.10457)

    def matProp_etangential_t1350(phi,baf):
        if(baf==baf_tan_array[0]): return matProp_eiso(phi,-1.42840,-0.19563,0.18991,-0.02591)
        elif(baf==baf_tan_array[1]): return matProp_eiso(phi,-2.24680,0.48243,-0.07687,0.00464)
        elif(baf==baf_tan_array[2]): return matProp_eiso(phi,-2.82930,0.76088,-0.22314,0.02431)
        elif(baf==baf_tan_array[3]): return matProp_eiso(phi,-3.25550,0.90423,-0.33175,0.04329)
        elif(baf==baf_tan_array[4]): return matProp_eiso(phi,-4.44780,1.60320,-0.58683,0.07458)
        elif(baf==baf_tan_array[5]): return matProp_eiso(phi,-5.67140,2.41920,-0.86155,0.10668)


    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    ax.set_xlim([0,4])

    ax.plot([0,4],[0,0],'-k',lw=0.75,c='dimgray')
    ax.plot(phi_array,matProp_eiso(phi_array,-1.24080,0.00175,0.08533,-0.01253),'-k',c=cjam[2],label=r"$T=600^\circ$C")
    ax.plot(phi_array,matProp_eiso(phi_array,-1.52390,0.13048,0.06299,-0.01072),'--k',c=cjam[1],label=r"$T=1032^\circ$C")
    ax.plot(phi_array,matProp_eiso(phi_array,-1.42840,-0.19563,0.18991,-0.02591),'-.k',c=cjam[0],label=r"$T=1350^\circ$C")

#    ax.xaxis.set_major_locator(MultipleLocator(300))
    ax.minorticks_on()


    ax.set_xlabel(r"$\phi$ (x10$^{25}$ n/m$^2$)")
    ax.set_ylabel(r"$\epsilon_{iso}(\rho_0)$ (%)")

    ax.legend(loc="lower left")

    plt.savefig('../figures/plot_pyc_mech_v4.pdf', bbox_inches='tight');
    plt.close(fig)

    fig = plt.figure(figsize=[11,5])
    gs = gridspec.GridSpec(1,2, left=0.2, right=1.3, top=1.2, bottom=0, wspace= 0.05)
    ax = fig.add_subplot(gs[0], projection='3d')

    editorial_settings.make_panes_transparent(ax)

    ax.set_xlim([0,4])

    X, Y = np.meshgrid(phi_array,baf_array);
    Z600 = []; Z1032 = []; Z1350 = [];
    for i in range(len(baf_array)):
        Z600.append([])
        Z1032.append([])
        Z1350.append([])
        [Z600[i].append( matProp_eradial_t600(phi_array[j],baf_array[i]) ) for j in range(len(phi_array))]
        [Z1032[i].append( matProp_eradial_t1032(phi_array[j],baf_array[i]) ) for j in range(len(phi_array))]
        [Z1350[i].append( matProp_eradial_t1350(phi_array[j],baf_array[i]) ) for j in range(len(phi_array))]
    Z600= np.array(Z600)
    Z1032= np.array(Z1032)
    Z1350= np.array(Z1350)

    ax.plot_surface(X, Y, Z600, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)
    ax.plot_surface(X, Y, Z1032, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)
    ax.plot_surface(X, Y, Z1350, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)

    ax.plot(phi_array, np.ones(len(phi_array))*baf_array[0], matProp_eradial_t600(phi_array,baf_array[0]), ':k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*baf_array[len(baf_array)-1], matProp_eradial_t600(phi_array,baf_array[len(baf_array)-1]) ,':k', lw=1)

    ax.plot(phi_array, np.ones(len(phi_array))*baf_array[0], matProp_eradial_t1032(phi_array,baf_array[0]) ,'--k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*baf_array[len(baf_array)-1], matProp_eradial_t1032(phi_array,baf_array[len(baf_array)-1]) ,'--k', lw=1)

    ax.plot(phi_array, np.ones(len(phi_array))*baf_array[0], matProp_eradial_t1350(phi_array,baf_array[0]) ,'-k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*baf_array[len(baf_array)-1], matProp_eradial_t1350(phi_array,baf_array[len(baf_array)-1]) ,'-k', lw=1)

    ax.plot(np.ones(len(baf_array))*phi_array[0], baf_array, np.zeros(len(baf_array)),'-k.', lw=1, ms=3)

    zx600= []
    [zx600.append(matProp_eradial_t600(phi_array[len(phi_array)-1],baf_array[i])) for i in range(len(baf_array))]
    ax.plot(np.ones(len(baf_array))*phi_array[len(phi_array)-1], baf_array, zx600,':k.', lw=1, ms=3, label=r"T=600$^\circ$C")

    zx1032= []
    [zx1032.append(matProp_eradial_t1032(phi_array[len(phi_array)-1],baf_array[i])) for i in range(len(baf_array))]
    ax.plot(np.ones(len(baf_array))*phi_array[len(phi_array)-1], baf_array, zx1032,'--k.', lw=1, ms=3, label=r"T=1032$^\circ$C")

    zx1350= []
    [zx1350.append(matProp_eradial_t1350(phi_array[len(phi_array)-1],baf_array[i])) for i in range(len(baf_array))]
    ax.plot(np.ones(len(baf_array))*phi_array[len(phi_array)-1], baf_array, zx1350,'-k.', lw=1, ms=3, label=r"T=1350$^\circ$C")

    ax.set_xlabel(r"$\phi$ (x10$^{25}$ n/m$^2$)")
    ax.set_ylabel(r"$BAF$ (-)")
    ax.set_zlabel(r"$\epsilon_{radial}(\rho_0)$ (%)")

    ax.legend(bbox_to_anchor=(0.25, 0.0), ncol=3)

    plt.savefig('../figures/plot_pyc_mech_v5.pdf', bbox_inches='tight');
    plt.close(fig)




    fig = plt.figure(figsize=[11,5])
    gs = gridspec.GridSpec(1,2, left=0.2, right=1.3, top=1.2, bottom=0, wspace= 0.05)
    ax = fig.add_subplot(gs[0], projection='3d')

    editorial_settings.make_panes_transparent(ax)

    ax.set_xlim([0,4])

    X, Y = np.meshgrid(phi_array,baf_tan_array);
    Z600 = []; Z1032 = []; Z1350 = [];
    for i in range(len(baf_tan_array)):
        Z600.append([])
        Z1032.append([])
        Z1350.append([])
        [Z600[i].append( matProp_etangential_t600(phi_array[j],baf_tan_array[i]) ) for j in range(len(phi_array))]
        [Z1032[i].append( matProp_etangential_t1032(phi_array[j],baf_tan_array[i]) ) for j in range(len(phi_array))]
        [Z1350[i].append( matProp_etangential_t1350(phi_array[j],baf_tan_array[i]) ) for j in range(len(phi_array))]
    Z600= np.array(Z600)
    Z1032= np.array(Z1032)
    Z1350= np.array(Z1350)

    ax.plot_surface(X, Y, Z600, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)
    ax.plot_surface(X, Y, Z1032, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)
    ax.plot_surface(X, Y, Z1350, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)

    ax.plot(phi_array, np.ones(len(phi_array))*baf_tan_array[0], matProp_etangential_t600(phi_array,baf_tan_array[0]), ':k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*baf_tan_array[len(baf_tan_array)-1], matProp_etangential_t600(phi_array,baf_tan_array[len(baf_tan_array)-1]) ,':k', lw=1)

    ax.plot(phi_array, np.ones(len(phi_array))*baf_tan_array[0], matProp_etangential_t1032(phi_array,baf_tan_array[0]) ,'--k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*baf_tan_array[len(baf_tan_array)-1], matProp_etangential_t1032(phi_array,baf_tan_array[len(baf_tan_array)-1]) ,'--k', lw=1)

    ax.plot(phi_array, np.ones(len(phi_array))*baf_tan_array[0], matProp_etangential_t1350(phi_array,baf_tan_array[0]) ,'-k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*baf_tan_array[len(baf_tan_array)-1], matProp_etangential_t1350(phi_array,baf_tan_array[len(baf_tan_array)-1]) ,'-k', lw=1)

    ax.plot(np.ones(len(baf_tan_array))*phi_array[0], baf_tan_array, np.zeros(len(baf_tan_array)),'-k.', lw=1, ms=3)

    zx600= []
    [zx600.append(matProp_etangential_t600(phi_array[len(phi_array)-1],baf_tan_array[i])) for i in range(len(baf_tan_array))]
    ax.plot(np.ones(len(baf_tan_array))*phi_array[len(phi_array)-1], baf_tan_array, zx600,':k.', lw=1, ms=3, label=r"600$^\circ$C")

    zx1032= []
    [zx1032.append(matProp_etangential_t1032(phi_array[len(phi_array)-1],baf_tan_array[i])) for i in range(len(baf_tan_array))]
    ax.plot(np.ones(len(baf_tan_array))*phi_array[len(phi_array)-1], baf_tan_array, zx1032,'--k.', lw=1, ms=3, label=r"1032$^\circ$C")

    zx1350= []
    [zx1350.append(matProp_etangential_t1350(phi_array[len(phi_array)-1],baf_tan_array[i])) for i in range(len(baf_tan_array))]
    ax.plot(np.ones(len(baf_tan_array))*phi_array[len(phi_array)-1], baf_tan_array, zx1350,'-k.', lw=1, ms=3, label=r"1350$^\circ$C")

    ax.set_xlabel(r"$\phi$ (x10$^{25}$ n/m$^2$)")
    ax.set_ylabel(r"$BAF$ (-)")
    ax.set_zlabel(r"$\epsilon_{tangential}(\rho_0)$ (%)")

    ax.legend(bbox_to_anchor=(0.25, 0.0), ncol=3)

    plt.savefig('../figures/plot_pyc_mech_v6.pdf', bbox_inches='tight');
    plt.close(fig)


    def matProp_alphar(tk,baf): return (30-37.5*2/(2+baf)) * (1+0.11*(tk-673)/700)
    def matProp_alphat(tk,baf): return (36*np.power(((1+baf)/(2+baf)-1),2)+1) * (1+0.11*(tk-673)/700)

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    ax.set_xlim([300,3000])
    ax.xaxis.set_major_locator(MultipleLocator(300))

    ax.plot([0,3000],[0,0],'-k',lw=0.75,c='dimgray')

    bafx= [1.1, 1.2, 1.3]
    for i in range(len(bafx)):
        ax.plot(tk_array,matProp_alphar(tk_array,bafx[i]),'-.k',c=cjam[i+1], dashes=[4,2],label=r"radial BAF="+str(bafx[i]))
    for i in range(len(bafx)):
        ax.plot(tk_array,matProp_alphat(tk_array,bafx[i]),'--k',c=cjam[i+1], dashes=[8,2],label=r"tangential BAF="+str(bafx[i]))

    ax.plot(tk_array,matProp_alphar(tk_array,1),'-k')

#    ax.xaxis.set_major_locator(MultipleLocator(300))
    ax.minorticks_on()

    ax.text(2400,6.5,r"BAF=1.0",rotation=9)


    ax.set_xlabel(r"$T$ (K)")
    ax.set_ylabel(r"$\alpha$ (10$^{-6}$/K)")

    ax.legend(loc="lower left",ncol=2)

    plt.savefig('../figures/plot_pyc_mech_v3.pdf', bbox_inches='tight');
    plt.close(fig)


    def matProp_k0(tk): return 2.193e-29-4.85e-32*(tk-273.15)+4.0147e-35*np.power((tk-273.15),2)
    def matProp_kss(tk,rho): return matProp_k0(tk) * (1+2.38*(1.9-rho)) * 2.0

    fig = plt.figure(figsize=[6.5,5.5])
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    ax.set_xlim([300,3000])
    ax.xaxis.set_major_locator(MultipleLocator(300))

#    ax.plot(tk_array,matProp_k0(tk_array),'-k',c=cjam[0],label=r"Eq.(3.28)")

    ax.plot(tk_array,matProp_kss(tk_array,1.7),'-k',c=cjam[2],label=r"$\rho=1.7$g/cm$^3$")
    ax.plot(tk_array,matProp_kss(tk_array,1.8),'-k',c=cjam[1],label=r"$\rho=1.8$g/cm$^3$")
    ax.plot(tk_array,matProp_kss(tk_array,1.9),'-k',c=cjam[0],label=r"$\rho=1.9$g/cm$^3$")
    ax.plot(tk_array,matProp_kss(tk_array,2.0),'-k',c=cjam[1],label=r"$\rho=2.0$g/cm$^3$")
    ax.plot(tk_array,matProp_kss(tk_array,2.1),'-k',c=cjam[2],label=r"$\rho=2.1$g/cm$^3$")

    ax.plot([0,3000],[0,0],'-k',lw=0.75,c='dimgray')

    ax.text(2400,1.125e-28,r"$\rho=2.1$g/cm$^3$",rotation=31)
    ax.text(2400,1.65e-28,r"$\rho=2.0$g/cm$^3$",rotation=42)
    ax.text(2400,2.2e-28,r"$\rho=1.9$g/cm$^3$",rotation=48)
    ax.text(2400,2.75e-28,r"$\rho=1.8$g/cm$^3$",rotation=52)
    ax.text(2400,3.3e-28,r"$\rho=1.7$g/cm$^3$",rotation=57)

    ax.minorticks_on()

    ax.set_xlabel(r"$T$ (K)")
    ax.set_ylabel(r"$K$ (-)")

    ax.legend(loc="upper left",ncol=1)

    plt.savefig('../figures/plot_pyc_mech_v2.pdf', bbox_inches='tight');
    plt.close(fig)

    def matProp_Er(tk,phi,rho=1.9,baf=1.0,Lc=30): return 25.5*(0.384+0.324e-3*rho)*(1.463-0.463*baf)*(2.985-0.0662*Lc)*(1+0.23*phi)*(1+0.00015*((tk-273.15)-20))
    def matProp_Et(tk,phi,rho=1.9,baf=1.0,Lc=30): return 25.5*(0.384+0.324e-3*rho)*(0.481+0.519*baf)*(2.985-0.0662*Lc)*(1+0.23*phi)*(1+0.00015*((tk-273.15)-20))
    def matProp_Etot(tk,phi,rho=1.9,baf=1.0,Lc=30): return 0.5 * ( matProp_Er(tk,phi,rho,baf,Lc) + matProp_Et(tk,phi,rho,baf,Lc) )

    fig = plt.figure(figsize=[11,5])
    gs = gridspec.GridSpec(1,2, left=0.2, right=1.3, top=1.2, bottom=0, wspace= 0.05)
    ax = fig.add_subplot(gs[0], projection='3d')

    editorial_settings.make_panes_transparent(ax)

    ax.set_xlim([0,4])

    X, Y = np.meshgrid(phi_array,tk_array);
    Z600 = []; Z1032 = []; Z1350 = [];
    for i in range(len(tk_array)):
        Z600.append([])
        Z1032.append([])
        Z1350.append([])
        [Z600[i].append( matProp_Er(tk_array[i],phi_array[j],baf=1.0) ) for j in range(len(phi_array))]
        [Z1032[i].append( matProp_Er(tk_array[i],phi_array[j],baf=1.1) ) for j in range(len(phi_array))]
        [Z1350[i].append( matProp_Er(tk_array[i],phi_array[j],baf=1.2) ) for j in range(len(phi_array))]
    Z600= np.array(Z600)
    Z1032= np.array(Z1032)
    Z1350= np.array(Z1350)
#
    ax.plot_surface(X, Y, Z600, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)
    ax.plot_surface(X, Y, Z1032, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)
    ax.plot_surface(X, Y, Z1350, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)
#
    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[0], matProp_Er(tk_array[0],phi_array,baf=1.0), ':k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[len(tk_array)-1], matProp_Er(tk_array[len(tk_array)-1],phi_array,baf=1.0), ':k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[0], tk_array, matProp_Er(tk_array,phi_array[0],baf=1.0), ':k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[len(phi_array)-1], tk_array, matProp_Er(tk_array,phi_array[len(phi_array)-1],baf=1.0), ':k', lw=1, label=r"BAF=1.0")

    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[0], matProp_Er(tk_array[0],phi_array,baf=1.1), '--k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[len(tk_array)-1], matProp_Er(tk_array[len(tk_array)-1],phi_array,baf=1.1), '--k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[0], tk_array, matProp_Er(tk_array,phi_array[0],baf=1.1), '--k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[len(phi_array)-1], tk_array, matProp_Er(tk_array,phi_array[len(phi_array)-1],baf=1.1), '--k', lw=1, label=r"BAF=1.1")

    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[0], matProp_Er(tk_array[0],phi_array,baf=1.2), '-k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[len(tk_array)-1], matProp_Er(tk_array[len(tk_array)-1],phi_array,baf=1.2), '-k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[0], tk_array, matProp_Er(tk_array,phi_array[0],baf=1.2), '-k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[len(phi_array)-1], tk_array, matProp_Er(tk_array,phi_array[len(phi_array)-1],baf=1.2), '-k', lw=1, label=r"BAF=1.2")

    ax.set_xlabel(r"$\phi$ (x10$^{25}$ n/m$^2$)")
    ax.set_ylabel(r"$T$ (K)")
    ax.set_zlabel(r"$E_r$ (GPa)")

    ax.yaxis.set_major_locator(MultipleLocator(300))
    ax.zaxis.set_major_locator(MultipleLocator(2))

    ax.legend(bbox_to_anchor=(0.25, 0.0), ncol=3)

    ax.set_zlim([10,30])

    plt.savefig('../figures/plot_pyc_mech_v1a.pdf', bbox_inches='tight');
    plt.close(fig)




    fig = plt.figure(figsize=[11,5])
    gs = gridspec.GridSpec(1,2, left=0.2, right=1.3, top=1.2, bottom=0, wspace= 0.05)
    ax = fig.add_subplot(gs[0], projection='3d')

    editorial_settings.make_panes_transparent(ax)

    ax.set_xlim([0,4])

    X, Y = np.meshgrid(phi_array,tk_array);
    Z600 = []; Z1032 = []; Z1350 = [];
    for i in range(len(tk_array)):
        Z600.append([])
        Z1032.append([])
        Z1350.append([])
        [Z600[i].append( matProp_Et(tk_array[i],phi_array[j],baf=1.0) ) for j in range(len(phi_array))]
        [Z1032[i].append( matProp_Et(tk_array[i],phi_array[j],baf=1.1) ) for j in range(len(phi_array))]
        [Z1350[i].append( matProp_Et(tk_array[i],phi_array[j],baf=1.2) ) for j in range(len(phi_array))]
    Z600= np.array(Z600)
    Z1032= np.array(Z1032)
    Z1350= np.array(Z1350)
#
    ax.plot_surface(X, Y, Z600, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)
    ax.plot_surface(X, Y, Z1032, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)
    ax.plot_surface(X, Y, Z1350, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)
#
    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[0], matProp_Et(tk_array[0],phi_array,baf=1.0), ':k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[len(tk_array)-1], matProp_Et(tk_array[len(tk_array)-1],phi_array,baf=1.0), ':k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[0], tk_array, matProp_Et(tk_array,phi_array[0],baf=1.0), ':k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[len(phi_array)-1], tk_array, matProp_Et(tk_array,phi_array[len(phi_array)-1],baf=1.0), ':k', lw=1, label=r"BAF=1.0")

    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[0], matProp_Et(tk_array[0],phi_array,baf=1.1), '--k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[len(tk_array)-1], matProp_Et(tk_array[len(tk_array)-1],phi_array,baf=1.1), '--k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[0], tk_array, matProp_Et(tk_array,phi_array[0],baf=1.1), '--k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[len(phi_array)-1], tk_array, matProp_Et(tk_array,phi_array[len(phi_array)-1],baf=1.1), '--k', lw=1, label=r"BAF=1.1")

    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[0], matProp_Et(tk_array[0],phi_array,baf=1.2), '-k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[len(tk_array)-1], matProp_Et(tk_array[len(tk_array)-1],phi_array,baf=1.2), '-k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[0], tk_array, matProp_Et(tk_array,phi_array[0],baf=1.2), '-k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[len(phi_array)-1], tk_array, matProp_Et(tk_array,phi_array[len(phi_array)-1],baf=1.2), '-k', lw=1, label=r"BAF=1.2")

    ax.set_xlabel(r"$\phi$ (x10$^{25}$ n/m$^2$)")
    ax.set_ylabel(r"$T$ (K)")
    ax.set_zlabel(r"$E_t$ (GPa)")

    ax.yaxis.set_major_locator(MultipleLocator(300))
    ax.zaxis.set_major_locator(MultipleLocator(2))

    ax.legend(bbox_to_anchor=(0.25, 0.0), ncol=3)

    ax.set_zlim([10,30])

    plt.savefig('../figures/plot_pyc_mech_v1b.pdf', bbox_inches='tight');
    plt.close(fig)



    fig = plt.figure(figsize=[11,5])
    gs = gridspec.GridSpec(1,2, left=0.2, right=1.3, top=1.2, bottom=0, wspace= 0.05)
    ax = fig.add_subplot(gs[0], projection='3d')

    editorial_settings.make_panes_transparent(ax)

    ax.set_xlim([0,4])

    X, Y = np.meshgrid(phi_array,tk_array);
    Z600 = []; Z1032 = []; Z1350 = [];
    for i in range(len(tk_array)):
        Z600.append([])
        Z1032.append([])
        Z1350.append([])
        [Z600[i].append( matProp_Etot(tk_array[i],phi_array[j],baf=1.0) ) for j in range(len(phi_array))]
        [Z1032[i].append( matProp_Etot(tk_array[i],phi_array[j],baf=1.1) ) for j in range(len(phi_array))]
        [Z1350[i].append( matProp_Etot(tk_array[i],phi_array[j],baf=1.2) ) for j in range(len(phi_array))]
    Z600= np.array(Z600)
    Z1032= np.array(Z1032)
    Z1350= np.array(Z1350)
#
    ax.plot_surface(X, Y, Z600, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)
    ax.plot_surface(X, Y, Z1032, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)
    ax.plot_surface(X, Y, Z1350, cmap='binary', alpha=0.2, linewidth=0.4, antialiased=True, zorder=0)
#
    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[0], matProp_Etot(tk_array[0],phi_array,baf=1.0), ':k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[len(tk_array)-1], matProp_Etot(tk_array[len(tk_array)-1],phi_array,baf=1.0), ':k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[0], tk_array, matProp_Etot(tk_array,phi_array[0],baf=1.0), ':k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[len(phi_array)-1], tk_array, matProp_Etot(tk_array,phi_array[len(phi_array)-1],baf=1.0), ':k', lw=1, label=r"BAF=1.0")

    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[0], matProp_Etot(tk_array[0],phi_array,baf=1.1), '--k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[len(tk_array)-1], matProp_Etot(tk_array[len(tk_array)-1],phi_array,baf=1.1), '--k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[0], tk_array, matProp_Etot(tk_array,phi_array[0],baf=1.1), '--k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[len(phi_array)-1], tk_array, matProp_Etot(tk_array,phi_array[len(phi_array)-1],baf=1.1), '--k', lw=1, label=r"BAF=1.1")

    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[0], matProp_Etot(tk_array[0],phi_array,baf=1.2), '-k', lw=1)
    ax.plot(phi_array, np.ones(len(phi_array))*tk_array[len(tk_array)-1], matProp_Etot(tk_array[len(tk_array)-1],phi_array,baf=1.2), '-k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[0], tk_array, matProp_Etot(tk_array,phi_array[0],baf=1.2), '-k', lw=1)
    ax.plot(np.ones(len(tk_array))*phi_array[len(phi_array)-1], tk_array, matProp_Etot(tk_array,phi_array[len(phi_array)-1],baf=1.2), '-k', lw=1, label=r"BAF=1.2")

    ax.set_xlabel(r"$\phi$ (x10$^{25}$ n/m$^2$)")
    ax.set_ylabel(r"$T$ (K)")
    ax.set_zlabel(r"$E_t$ (GPa)")

    ax.yaxis.set_major_locator(MultipleLocator(300))

    ax.legend(bbox_to_anchor=(0.175, 0.0), ncol=3)

    plt.savefig('../figures/plot_pyc_mech_v1c.pdf', bbox_inches='tight');
    plt.close(fig)
