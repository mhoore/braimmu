import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import axes3d
from math import sqrt, floor, pi
from matplotlib import rc, rcParams, cm
from matplotlib.ticker import MultipleLocator
import matplotlib.mlab as mlab
#from matplotlib.mlab import griddata
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import sys
from scipy.interpolate import griddata
import numpy.ma as ma
from matplotlib.colors import LogNorm
import subprocess
import os
import nibabel as nib

rc('font',**{'family':'sans-serif','serif':['Times'], 'weight': 'bold', 'size' : 24})
rc('text', usetex=True)
rcParams.update({'figure.autolayout': True})

command = "mkdir mean_field"
p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

# Data Input

if len(sys.argv) < 1:
    print(sys.argv)
    sys.exit("Error: chack the arguments")

fr = sys.argv[1]

dt = 0.001
delt = dt / 365.

T_id = 0
M_id = 1
N_id = 2
S_id = 3
F_id = 4
A_id = 5

names = (r'$\rm time ~{\rm (yr)}$', r'$M / M_{\rm max}$', r'$N / N_{\rm max}$', r'$S / S_{\rm max}$', r'$F / F_{\rm max}$', r'$\rm Astrogliosis$')
names_real = (r'$t ~{\rm (yr)}$', r'$M ~{\rm (mL^{-1})}$', r'$N ~{\rm (mL^{-1})}$', r'$S ~{\rm (\mu M)}$',r'$F ~{\rm (\mu M)}$',r'$\rm Astrogliosis$')

class Subplots:
    """ Arrange subplot grid structure (square box is assumed)"""

    totcnt = -1             # Total number of subplots

    # Constructor
    def __init__(self, f, lx, ly, sx, sy, bx, by, t):
        self.fig = f        # Figure axes handle
        self.length_x = lx  # Length of the subplot box
        self.length_y = ly  # Length of the subplot box
        self.sepx = sx      # Separation distance between subplots
        self.sepy = sy      # Separation distance between subplots

        self.begx = bx      # Beginning (offset) in the figure box
        self.begy = by      # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in x direction

    # Add a subplot in the grid structure
    def addSubplot(self):
        
        # Increase the number of subplots in the figure
        self.totcnt += 1
        
        # Indices of the subplot in the figure
        self.nx = self.totcnt%(self.tot)
        self.ny = int(self.totcnt/(self.tot))
        
        self.xbeg = self.begx + self.nx*self.length_x + self.nx*self.sepx
        self.ybeg = self.begy + self.ny*self.length_y + self.ny*self.sepy
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length_x,self.length_y])

def main():
    
    #data = np.genfromtxt(fr, dtype=np.float32, skip_header=2)
    # step region mic neu sAb fAb ast cir
    data = np.genfromtxt(fr, dtype=np.float32, skip_header=2, usecols=(0,2,3,4,5,6))
    
    data_PAR = data[0::2,:]
    data_CSF = data[1::2,:]
    
    #############
    ### plots
    #############
    fig = pl.figure(figsize=(7,5))
    ax = fig.add_subplot(111)
    
    # dementia (1 - N/N0)
    x = data_PAR[:,T_id] * delt
    y = 1.0 - data_PAR[:,N_id] / np.max(data_PAR[:,N_id])
    lbl = r'$\rm dementia$'
    
    ax.plot(x, y, c='k', ls='-', lw=2.0, label=lbl)

    # fAb
    x = data_PAR[:,T_id] * delt
    y = data_PAR[:,F_id] / np.max(data_PAR[:,F_id])
    lbl = r'$\rm F / F_{\rm max}$'
    
    ax.plot(x, y, c='r', ls='--', lw=2, label=lbl)

    # A
    x = data_PAR[:,T_id] * delt
    y = data_PAR[:,A_id]
    lbl = r'$\rm astrogliosis$'
    
    ax.plot(x, y, c='b', ls='-.', lw=2, label=lbl)

    # PAR sAb
    x = data_PAR[:,T_id] * delt
    y = data_PAR[:,S_id] / np.max(data_PAR[:,S_id])
    lbl = r'$\rm S / S_{\rm max} (ISF)$'
    
    ax.plot(x, y, c='g', ls='-', lw=2, label=lbl)

    # CSF sAb
    x = data_CSF[:,T_id] * delt
    y = data_CSF[:,S_id] / np.max(data_CSF[:,S_id])
    lbl = r'$\rm S / S_{\rm max} (CSF)$'
    
    ax.plot(x, y, c='g', ls='-.', lw=2, label=lbl)


    ax.legend(loc='upper center', borderpad=.4, labelspacing=.1, borderaxespad=.1, columnspacing=.2, fontsize=16, ncol=3)

    #ax.set_yscale('log')

    ax.set_xlim(0.0,np.max(x))
    ax.set_xticks([0,5,10])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel(names[T_id])
    ax.set_ylim(0,1.3)
    ax.set_yticks(np.linspace(0, 1, 3, endpoint=True))
    #ax.set_ylabel('y', color='black', fontsize = 24)

    minorLocator   = MultipleLocator(0.5)
    ax.xaxis.set_minor_locator(minorLocator)
    #minorLocator   = MultipleLocator(0.1)
    #ax.yaxis.set_minor_locator(minorLocator)

    pl.savefig("mean_field/mean_field.png", format = 'png', dpi=100, orientation='landscape')
    pl.show()
    #sys.stdin.read(1)
    pl.close()
    
if __name__ == '__main__':
    main()
