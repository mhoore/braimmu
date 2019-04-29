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

command = "mkdir single"
p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

# Data Input
scale_fac = 1.e-3

if len(sys.argv) < 7:
    print(sys.argv)
    sys.exit("Error: chack the arguments")

fr = sys.argv[1]
dim = int(sys.argv[2])
sec = int(sys.argv[3])
time = sys.argv[4]
realtime = float(sys.argv[5])
fr0 = sys.argv[6]

AGENTS = ("mic","neu","sAb","fAb","ast","typ")
ag_tit = (r'$\rm Microglia$', r'$\rm Neurons$', r'$\rm [sA\beta]$',r'$\rm [fA\beta]$',r'$\rm Astrogliosis$',r'$\rm Tissue$',r'$\rm Atrophy$')
ag_names = (r'$M ~{\rm (mL^{-1})}$', r'$N ~{\rm (mL^{-1})}$', r'$S ~{\rm (\mu M)}$',r'$F ~{\rm (\mu M)}$',r'$\rm Astrogliosis$',r'$\rm Tissue$',r'$\rm Atrophy$')
ag_sgn = ('M', 'N', 'S','F','A','TYPE','atrophy')

x = []
y = []
Npoints = []

agents = []
agents0 = []

Nx = np.zeros(4, dtype=np.int32)
dx = np.zeros(3, dtype=np.float32)
lx = np.zeros(3, dtype=np.float32)

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

    img0 = nib.load(fr0)
    data0 = img0.get_fdata()
    data_sec0 = []
    
    Nx[0], Nx[1], Nx[2], dum, Nx[3] = data.shape
    dx[0], dx[1], dx[2], dt, dum = img.header.get_zooms()
    
    for i in range(3):
        dx[i] *= scale_fac
        lx[i] = Nx[i] * dx[i]
    
    if (dim == 0):
        data_sec.append( data[sec-1,:,:,0,:] )
        data_sec0.append( data0[sec-1,:,:,0,:] )
        Npoints.append( Nx[1] * Nx[2] )
        x.append( np.zeros(Npoints[0], dtype=np.float32) )
        y.append( np.zeros(Npoints[0], dtype=np.float32) )

    if (dim == 1):
        data_sec.append( data[:,sec-1,:,0,:] )
        data_sec0.append( data0[:,sec-1,:,0,:] )
        Npoints.append( Nx[0] * Nx[2] )
        x.append( np.zeros(Npoints[0], dtype=np.float32) )
        y.append( np.zeros(Npoints[0], dtype=np.float32) )
    
    if (dim == 2):
        data_sec.append( data[:,:,sec-1,0,:] )
        data_sec0.append( data0[:,:,sec-1,0,:] )
        Npoints.append( Nx[0] * Nx[1] )
        x.append( np.zeros(Npoints[0], dtype=np.float32) )
        y.append( np.zeros(Npoints[0], dtype=np.float32) )

    for j in range(Nx[3]):
        agents.append( np.zeros(Npoints[0], dtype=np.float32) )
        agents0.append( np.zeros(Npoints[0], dtype=np.float32) )
    
    ii = 0
    jj = 0
    if (dim == 0):
        ii = 1
        jj = 2
    elif (dim == 1):
        ii = 0
        jj = 2
    elif (dim == 2):
        ii = 0
        jj = 1
    
    c = 0
    for i in range(Nx[ii]):
        dum = dx[ii] * i
        for j in range(Nx[jj]):
            x[0][c] = dum
            y[0][c] = dx[jj] * j
            
            for ag_id in range(Nx[3]):
                agents[ag_id][c] = data_sec[0][i,j,ag_id]
                agents0[ag_id][c] = data_sec0[0][i,j,ag_id]
            
            c += 1
    
    for me in range(Nx[3]-1):
        z = data[:,:,:,0,me].flatten()
        SMALL = 1.e-6
        z[z <= SMALL] = 0
        Z = np.ma.masked_equal(z,0)
        z = Z.compressed() # get normal array with masked values removed
        
        mu = np.mean(z)
        sig = np.std(z)
        
        print(me,ag_sgn[me],mu,sig)

        zs = (z-mu)/sig
        if (math.isnan(sig)):
            continue
        
        cut_min = np.min(zs)
        cut_max = np.max(zs)
        
        Zscore = ( agents[me] - mu)/sig
        
        fw = "single/s_Z_" + AGENTS[me] + "_d" + str(dim) + "_s" + str(sec) + "_t" + time + ".png"

        z   = data[:,:,:,0,me].flatten()
        tis = data[:,:,:,0,-1].flatten()

        make_plot(Zscore,z,tis,me,cut_min,cut_max,fw,100)

    # tissue
    me = Nx[3]-1
    z = data[:,:,:,0,me].flatten()
    Zscore = agents[me]
    fw = "single/s_type_d" + str(dim) + "_s" + str(sec) + "_t" + time + ".png"
    make_plot(Zscore,z,tis,me,-1,3,fw,4)

    # attrophy (1-N)
    me = 1
    atrophy = - agents[me] / agents0[me] + 1.0

    z = - data[:,:,:,0,me].flatten() / data0[:,:,:,0,me].flatten() + 1.0
    SMALL = 1.e-6
    z[z <= SMALL] = 0
    Z = np.ma.masked_equal(z,0)
    z = Z.compressed() # get normal array with masked values removed
    
    cut_min = 0.0
    cut_max = 1.0
    
    fw = "single/s_atrophy_d" + str(dim) + "_s" + str(sec) + "_t" + time + ".png"

    z = - data[:,:,:,0,me].flatten() / data0[:,:,:,0,me].flatten() + 1.0
    tis = data[:,:,:,0,-1].flatten()

    make_plot(atrophy,z,tis,Nx[3],cut_min,cut_max,fw,100)
    
    data = []
    data_sec = []
    data0 = []
    data_sec0 = []

## plots
def make_plot(Zscore,z,tis,me,cut_min,cut_max,fw,contour_levels):
    fig = pl.figure(figsize=(9,9))
    #ax = fig.add_subplot(111)
    ax = fig.add_axes([0.02,0.02,0.9,0.9])
    
    mlx   = MultipleLocator(10)
    mly   = MultipleLocator(10)

    lbl = r'$z_{%s}$' % ag_sgn[me]
    if (me == Nx[3]-1):
        lbl = r'$\rm Tissue$'
    elif (me == Nx[3]):
        lbl = r'$\rm Atrophy$'

    xls = []
    yls = []
    if (dim == 0):
        xls = np.linspace(0,lx[1],Nx[1])
        yls = np.linspace(0,lx[2],Nx[2])
    elif (dim == 1):
        xls = np.linspace(0,lx[0],Nx[0])
        yls = np.linspace(0,lx[2],Nx[2])
    elif (dim == 2):
        xls = np.linspace(0,lx[0],Nx[0])
        yls = np.linspace(0,lx[1],Nx[1])
    
    x_gr, y_gr = np.meshgrid(xls, yls)
        
    z_gr = griddata((x[0],y[0]), Zscore, (x_gr, y_gr), method='nearest')
    
    cs = ax.contourf(x_gr,y_gr,z_gr, contour_levels, cmap=pl.cm.jet,
                            levels = np.linspace(cut_min,cut_max,contour_levels), extend='both') # , norm = LogNorm())
    #cs = ax.pcolor(x_gr,y_gr,z_gr, cmap=pl.cm.jet, vmin=cut_min, vmax=cut_max)
    
    cs.cmap.set_under('None')
    cs.cmap.set_over('k')

    # scale bar
    from matplotlib.patches import Rectangle
    #cax = pl.gca()
    ax.add_patch(Rectangle((0, 0), 100, 2, alpha=1, color='k'))
    #fig.text(0.1,0.1,r'$\rm 100 ~mm$')

    # ax.plot(x, y, 'ko', markersize=4)

    if (dim == 0):
        ax.set_xlim(0,lx[1] )
        ax.set_ylim(0,lx[2] )
    elif (dim == 1):
        ax.set_xlim(0,lx[0] )
        ax.set_ylim(0,lx[2] )
    elif (dim == 2):
        ax.set_xlim(0,lx[0] )
        ax.set_ylim(0,lx[1] )

    ax.xaxis.set_minor_locator(mlx)
    ax.yaxis.set_minor_locator(mly)

    tit = r'%s ,  $\rm t = %g ~(day)$' % (ag_tit[me],realtime)
    ax.set_xticks([])
    ax.set_yticks([])

    fig.text(0.50,0.96,tit,fontsize=30, ha='center', va='center')
    
    cbar_ax = fig.add_axes([0.93,0.15,0.02,0.7])
    cbar = pl.colorbar(cs,cax=cbar_ax)
    cbar.set_ticks([cut_min,cut_max])
    cbar.ax.set_yticklabels([r'$\rm %.1f$' % cut_min, r'$\rm %.1f$' % cut_max],rotation=90,fontsize=30)
    fig.text(0.96,0.5,lbl, rotation='vertical', fontsize=30)
    
    ax.set_aspect('equal')

    pl.savefig(fw, format = 'png', dpi=100, orientation='landscape')
    pl.close()

if __name__ == '__main__':
    main()
