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
import math

rc('font',**{'family':'sans-serif','serif':['Times'], 'weight': 'bold', 'size' : 26})
rc('text', usetex=True)
rcParams.update({'figure.autolayout': True})

if len(sys.argv) < 6:
    print(sys.argv)
    sys.exit("Error: check the arguments")

sec = []
fr = sys.argv[1]
fr0 = sys.argv[2]
sec.append(int(sys.argv[3]))
sec.append(int(sys.argv[4]))
sec.append(int(sys.argv[5]))

x = []
y = []
Npoints = []

Nx = np.zeros(4, dtype=np.int32)

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
    
    #data = np.genfromtxt(fr, dtype=np.float32, skip_header=2, max_rows = 1)
    img = nib.load(fr)
    data = img.get_fdata()
    data_sec = []
    
    descrip = img.header['descrip'].astype(str).item()
    agents = descrip.split()

    img0 = nib.load(fr0)
    data0 = img0.get_fdata()
    data_sec0 = []

    descrip = img0.header['descrip'].astype(str).item()
    agents0 = descrip.split()
    
    if (data.shape != data0.shape):
        print("Nifti files shape mismatch!")
        exit(1)
    
    if (agents != agents0):
        print("Nifti files description mismatch!")
        exit(1)
    
    Nx[0], Nx[1], Nx[2], dum, Nx[3] = data.shape
    
    # dimensions
    diff = []
    diff.append(data[sec[0]-1,:,:,0,:] - data0[sec[0]-1,:,:,0,:])
    diff.append(data[:,sec[1]-1,:,0,:] - data0[:,sec[1]-1,:,0,:])
    diff.append(data[:,:,sec[2]-1,0,:] - data0[:,:,sec[2]-1,0,:])

    for i in range(3):
        for me in range(Nx[3]-1):
            print(i, agents[me], np.sum(diff[i][:,:,me]))

if __name__ == '__main__':
    main()
