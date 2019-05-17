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

command = "mkdir slice_hists"
p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

command = "mkdir slice_contours"
p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

# Data Input
scale_fac = 1.e-3

if len(sys.argv) < 5:
    print(sys.argv)
    sys.exit("Error: chack the arguments")

fr = sys.argv[1]
dim = int(sys.argv[2])
sec_id = int(sys.argv[3])
time = sys.argv[4]
realtime = sys.argv[5]

AGENTS = ("mic","neu","sAb","fAb","ast","typ")
ag_names = (r'$M ~{\rm (mL^{-1})}$', r'$N ~{\rm (mL^{-1})}$', r'$S ~{\rm (\mu M)}$',r'$F ~{\rm (\mu M)}$',r'$\rm Astrogliosis$',r'$\rm Tissue$')
ag_sgn = ('M', 'N', 'S','F','A','TYPE')

x = []
agents = []
Zscore = []
Pearson = []

xlo = []
xhi = []

xlim = np.zeros(2, dtype=np.float32)
ylim = np.zeros(2, dtype=np.float32)

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
    
    Nx[0], Nx[1], Nx[2], dum, Nx[3] = data.shape
    dx[0], dx[1], dx[2], dt, dum = img.header.get_zooms()
    
    for i in range(3):
        lx[i] = Nx[i] * dx[i] * scale_fac
    
    if (dim == 0):
        xlbl = r'$y ~{\rm (mm)}$'
        ylbl = r'$z ~{\rm (mm)}$'
        ii = 1
        jj = 2
        data_sec.append(data[sec_id-1,:,:,0,:])
        
    elif (dim == 1):
        xlbl = r'$x ~{\rm (mm)}$'
        ylbl = r'$z ~{\rm (mm)}$'
        ii = 0
        jj = 2
        data_sec.append(data[:,sec_id-1,:,0,:])

    elif (dim == 2):
        xlbl = r'$x ~{\rm (mm)}$'
        ylbl = r'$y ~{\rm (mm)}$'
        ii = 0
        jj = 1
        data_sec.append(data[:,:,sec_id-1,0,:])

    Npoints = Nx[ii] * Nx[jj]
    
    x.append(np.zeros(Npoints, dtype=np.float32))
    x.append(np.zeros(Npoints, dtype=np.float32))

    for i in range(Nx[3]):
        agents.append(np.zeros(Npoints, dtype=np.float32))
    
    c = 0
    for i in range(Nx[ii]):
        dum = dx[ii] * i
        for j in range(Nx[jj]):
            x[0][c] = dum
            x[1][c] = dx[jj] * j
            
            for ag_id in range(Nx[3]):
                agents[ag_id][c] = data_sec[0][i,j,ag_id]
            
            c += 1
    
    for i in range(2):
        x[i] *= scale_fac
        xlo.append(np.min(x[i]))
        xhi.append(np.max(x[i]))
        
    xlim[0] = xlo[0]
    xlim[1] = xhi[0]
    
    ylim[0] = xlo[1]
    ylim[1] = xhi[1]
    
    data = []
    data_sec = []
    
    prune_data()
    
    print(x[0].shape,x[1].shape,agents[0].shape)
    for i in range(len(x[0])):
        if (agents[1][i] <= 0.0):
            print(agents[1][i], agents[-1][i])

    for i in range(Nx[3]):
        Zscore.append( ( agents[i] - np.mean(agents[i]) ) / np.std(agents[i]) )
        #if (np.std(agents[i]) == 0):
        #   Zscore[i] = 0.0
        #  print(agents[i],Zscore[i])
    

    for me in range(Nx[3]-1):

        labels = (xlbl,ylbl,ag_names[me])
            
        fw = "slice_hists/s_hist_" + AGENTS[me] + "_d" + str(dim) + "_s" + str(sec_id) + "_t" + time + ".png"
        cut_min, cut_max = make_histogram(agents[me],labels,fw)
        
        fw = "slice_contours/s_cont_" + AGENTS[me] + "_d" + str(dim) + "_s" + str(sec_id) + "_t" + time + ".png"
        #make_contour(x[0],x[1],agents[me],Nx[ii],Nx[jj],lx[ii],lx[jj],labels,cut_min,cut_max,fw)
            
        # Z-scores
        lbl = r'$z_{%s}$' % ag_sgn[me]
        labels = (xlbl,ylbl,lbl)

        fw = "slice_contours/s_cont_Z_" + AGENTS[me] + "_d" + str(dim) + "_s" + str(sec_id) + "_t" + time + ".png"
        
        cut_max = np.max(Zscore[me]) #max( abs(np.min(Zscore[me])),np.max(Zscore[me]) )
        cut_min = np.min(Zscore[me]) #- cut_max
        make_contour(x[0],x[1],Zscore[me],Nx[ii],Nx[jj],lx[ii],lx[jj],labels,cut_min, cut_max,fw)
            
        # Pearson's r correlation
        for me2 in range(me+1,Nx[3]-1):

            if (me != 1):
                continue
            if (me2 !=3  and me2 != 4):
                continue

            Pearson = Zscore[me] * Zscore[me2] / (len(Zscore[me]))
            
            lbl = r'${\rm corr}(%s,%s)$' % (ag_sgn[me],ag_sgn[me2])
            labels = (xlbl,ylbl,lbl)
            
            fw = "slice_contours/s_cont_corr_" + AGENTS[me] + "_" + AGENTS[me2] + "_d" + str(dim) + "_s" + str(sec_id) + "_t" + time + ".png"
            
            cut_max = np.max(Pearson) # max( abs(np.min(Pearson)),np.max(Pearson) )
            cut_min = np.min(Pearson) # - cut_max

            make_contour(x[0],x[1],Pearson,Nx[ii],Nx[jj],lx[ii],lx[jj],labels,cut_min, cut_max,fw)

def prune_data():
    # find useless data
    marray = np.ma.masked_where(agents[-1] <= 0,agents[-1]).mask
    
    x[0] = np.ma.array(x[0], mask = marray).compressed()
    x[1] = np.ma.array(x[1], mask = marray).compressed()
    
    for i in range(Nx[3]):
        agents[i] = np.ma.array(agents[i], mask = marray).compressed()

def make_histogram(z,labels,fw):
    # find useless data
    SMALL = 1.e-6
    z[z <= SMALL] = 0
    Z = np.ma.masked_equal(z,0)
    z = Z.compressed() # get normal array with masked values removed
    
    cut_min = np.min(z)
    cut_max = np.max(z)
    
    #print(ag_names[me],np.min(z),np.max(z))
    
    fig = pl.figure(figsize=(7,5))
    ax = fig.add_subplot(111)
    ax.hist(z, histtype='stepfilled',bins=50,density=True, fc='#CCCCCC', alpha=0.5, edgecolor='black', linewidth=1.2)  # bins='auto'
    ax.set_xlabel(labels[2])
    ax.set_ylabel(r'$\rm probability$')
    
    #ax.set_xlim(cut_min,cut_max)
    #ax.set_xticklabels([r'%.2e' % cut_min, r'%.2e' % cut_max])
    
    pl.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    pl.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    pl.savefig(fw, format = 'png', dpi=100, orientation='landscape')
    pl.close()
    
    return cut_min, cut_max
    
## contour plots
def make_contour(x,y,z,Nx,Ny,lx,ly,labels,cut_min,cut_max,fw):
    contour_levels = 100
    # Plot properties
    if (dim == 0):
        figsizex = 10.0
        figsizey =  8.0
    elif (dim == 1):
        figsizex = 10.0
        figsizey =  9.5
    elif (dim == 2):
        figsizex = 10.0
        figsizey = 11.0
    
    ar = figsizex / figsizey
    n_X = 1
    n_Y = 1
    ax_len_x = 0.78/n_X             # Length of one subplot square box
    ax_len_y = ax_len_x*ly/lx * ar  # Length of one subplot square box
    ax_bx = 0.10                    # Beginning/offset of the subplot in the box
    ax_by = 0.09 * ar               # Beginning/offset of the subplot in the box
    ax_sepx = 0.01                  # Separation length between two subplots
    ax_sepy = 0.02 * ar             # Separation length between two subplots
    total_subplots_in_x = n_X       # Total number of subplots in x direction

    ax = []
    fig = pl.figure(figsize=(figsizex,figsizey))
    subp = Subplots(fig, ax_len_x, ax_len_y, ax_sepx, ax_sepy, ax_bx, ax_by, total_subplots_in_x)
    subcnt = 0

    mlx   = MultipleLocator(10)
    mly   = MultipleLocator(10)

    xls = np.linspace(np.min(x),np.max(x),Nx)
    yls = np.linspace(np.min(y),np.max(y),Ny)
    x_gr, y_gr = np.meshgrid(xls, yls)
    
    for i_y in range(n_Y) :
        for i_x in range(n_X) :
            ax.append( subp.addSubplot() )
            
            # grid the data.
            z_gr = griddata((x,y), z, (x_gr, y_gr), method='cubic')

            #cs = ax[subcnt].contourf(x_gr,y_gr,z_gr, contour_levels, cmap=pl.cm.jet,
             #                        levels = np.linspace(cut_min,cut_max,contour_levels), extend='both') # , norm = LogNorm())
            cs = ax[subcnt].pcolor(x_gr,y_gr,z_gr, cmap=pl.cm.jet, vmin=cut_min, vmax=cut_max)
            
            cs.cmap.set_under('None')
            cs.cmap.set_over('k')

            # ax[subcnt].plot(x, y, 'ko', markersize=4)

            ax[subcnt].set_xlim(xlim[0],xlim[1])
            ax[subcnt].set_ylim(ylim[0],ylim[1])

            ax[subcnt].tick_params(axis='x',which='minor',bottom='on')
            ax[subcnt].tick_params(axis='y',which='minor',bottom='on')

            ax[subcnt].xaxis.set_minor_locator(mlx)
            ax[subcnt].yaxis.set_minor_locator(mly)

            tit = r'$\rm time = %s ~(day)$' % realtime
            #ax[subcnt].set_xticks([pi/2,pi,3.0*pi/2,2.0*pi])
            #ax[subcnt].set_yticks([pi/2,pi])

            #if (i_y == 1):
            ax[subcnt].set_title(tit)
            #ax[subcnt].set_xticklabels([])
            #if (i_y == 0):
            #ax[subcnt].set_xlabel(r'$\left( \phi - \pi \right) \sin{\theta} + \pi$')
            #ax[subcnt].set_xticklabels([r'$\rm \pi/2$',r'$\rm \pi$',r'$\rm 3\pi/2$',r'$\rm 2 \pi$'])

            ax[subcnt].set_xlabel(labels[0])
            ax[subcnt].set_ylabel(labels[1])
            
            #ax[subcnt].set_yticklabels([r'$\rm \pi/2$',r'$\rm \pi$'])
            #else:
            #pl.setp(ax[subcnt].get_yticklabels(), visible=False)
            #ax[subcnt].set_yticklabels([])

            #if (i_x == n_X - 1 and i_y == 0) :
            cbar_ax = fig.add_axes([n_X*ax_len_x + (n_X-1.0)*ax_sepx + ax_bx + 0.01, 
                                    ax_by + i_y*ax_sepy + (0.1 + i_y)*ax_len_y, 0.03, 0.8*ax_len_y])
            cbar = pl.colorbar(cs,cax=cbar_ax)
            cbar.set_ticks([cut_min,cut_max])
            #mticks = cbar.norm([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20])
            #cbar.ax.yaxis.set_ticks(mticks, minor=True)
            #cbar.set_ticklabels([0, r'%.2e' % np.max(z)])
            cbar.ax.set_yticklabels([r'%.2e' % cut_min, r'%.2e' % cut_max],rotation=90)
            fig.text(n_X*ax_len_x + (n_X-1.0)*ax_sepx + ax_bx + 0.06, 
                     ax_by + i_y*ax_sepy + (0.5 + i_y)*ax_len_y, labels[2], rotation='vertical')

            subcnt += 1

    pl.savefig(fw, format = 'png', dpi=100, orientation='landscape')
    pl.close()

if __name__ == '__main__':
    main()
