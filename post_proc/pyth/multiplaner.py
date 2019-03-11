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

rc('font',**{'family':'sans-serif','serif':['Times'], 'weight': 'bold', 'size' : 30})
rc('text', usetex=True)
rcParams.update({'figure.autolayout': True})

command = "mkdir multi"
p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

# Data Input
scale_fac = 1.e-3

if len(sys.argv) < 7:
    print(sys.argv)
    sys.exit("Error: chack the arguments")

sec = []
fr = sys.argv[1]
sec.append(int(sys.argv[2]))
sec.append(int(sys.argv[3]))
sec.append(int(sys.argv[4]))
time = sys.argv[5]
realtime = sys.argv[6]

AGENTS = ("mic","neu","sAb","fAb","ast","typ")
ag_names = (r'$M ~{\rm (mL^{-1})}$', r'$N ~{\rm (mL^{-1})}$', r'$S ~{\rm (\mu M)}$',r'$F ~{\rm (\mu M)}$',r'$\rm Astrogliosis$',r'$\rm Tissue$')
ag_sgn = ('M', 'N', 'S','F','A','TYPE')

x = []
y = []
Npoints = []

agents = []

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
        dx[i] *= scale_fac
        lx[i] = Nx[i] * dx[i]
    
    # dim0
    data_sec.append(data[sec[0]-1,:,:,0,:])
    Npoints.append(Nx[1] * Nx[2])
    x.append(np.zeros(Npoints[0], dtype=np.float32))
    y.append(np.zeros(Npoints[0], dtype=np.float32))
    agents.append([])
    
    # dim1
    data_sec.append(data[:,sec[1]-1,:,0,:])
    Npoints.append(Nx[0] * Nx[2])
    x.append(np.zeros(Npoints[1], dtype=np.float32))
    y.append(np.zeros(Npoints[1], dtype=np.float32))
    agents.append([])
    
    # dim2
    data_sec.append(data[:,:,sec[2]-1,0,:])
    Npoints.append(Nx[0] * Nx[1])
    x.append(np.zeros(Npoints[2], dtype=np.float32))
    y.append(np.zeros(Npoints[2], dtype=np.float32))
    agents.append([])
    
    for i in range(3):
        for j in range(Nx[3]):
            agents[i].append(np.zeros(Npoints[i], dtype=np.float32))
    
    # dim 0
    ii = 1
    jj = 2
    c = 0
    for i in range(Nx[ii]):
        dum = dx[ii] * i
        for j in range(Nx[jj]):
            x[0][c] = dum
            y[0][c] = dx[jj] * j
            
            for ag_id in range(Nx[3]):
                agents[0][ag_id][c] = data_sec[0][i,j,ag_id]
            
            c += 1
    
    # dim 1
    ii = 0
    jj = 2
    c = 0
    for i in range(Nx[ii]):
        dum = dx[ii] * i
        for j in range(Nx[jj]):
            x[1][c] = dum
            y[1][c] = dx[jj] * j
            
            for ag_id in range(Nx[3]):
                agents[1][ag_id][c] = data_sec[1][i,j,ag_id]
            
            c += 1
    
    # dim 2
    ii = 0
    jj = 1
    c = 0
    for i in range(Nx[ii]):
        dum = dx[ii] * i
        for j in range(Nx[jj]):
            x[2][c] = dum
            y[2][c] = dx[jj] * j
            
            for ag_id in range(Nx[3]):
                agents[2][ag_id][c] = data_sec[2][i,j,ag_id]
            
            c += 1
    
    prune_data()
    
    for me in range(Nx[3]-1):
        z = data[:,:,:,0,me].flatten()
        SMALL = 1.e-6
        z[z <= SMALL] = 0
        Z = np.ma.masked_equal(z,0)
        z = Z.compressed() # get normal array with masked values removed
        
        mu = np.mean(z)
        sig = np.std(z)
        
        print(mu,sig)

        zs = (z-mu)/sig
        
        cut_min = np.min(zs)
        cut_max = np.max(zs)
        
        Zscore = []
        for i in range(3):
            Zscore.append( ( agents[i][me] - mu)/sig )
        
        fw = "multi/m_Z_" + AGENTS[me] + "_s" + str(sec[0]) + "_" + str(sec[1]) + "_" + str(sec[2]) + "_t" + time + ".png"
        make_plot(Zscore,z,me,cut_min,cut_max,fw)

    data = []
    data_sec = []

def prune_data():
    # find useless data
    for i in range(3):
        marray = np.ma.masked_where(agents[i][-1] <= 0,agents[i][-1]).mask
        
        x[i] = np.ma.array(x[i], mask = marray).compressed()
        y[i] = np.ma.array(y[i], mask = marray).compressed()
    
        for j in range(Nx[3]):
            agents[i][j] = np.ma.array(agents[i][j], mask = marray).compressed()
    
## plots
def make_plot(Zscore,z,me,cut_min,cut_max,fw):
    contour_levels = 100
    fig = pl.figure(figsize=(15,14))
    ax = fig.add_subplot(111)
    
    xsep = 1.
    ysep = 1.
    
    mlx   = MultipleLocator(10)
    mly   = MultipleLocator(10)

    lbl = r'$z_{%s}$' % ag_sgn[me]

    ## dim0
    bx = lx[0] + xsep
    by = lx[1] + ysep
    xls = np.linspace(bx,bx+lx[1],Nx[1])
    yls = np.linspace(by,by+lx[2],Nx[2])
    x_gr, y_gr = np.meshgrid(xls, yls)
    
    # grid the data.
    z_gr = griddata((x[0]+bx,y[0]+by), Zscore[0], (x_gr, y_gr), method='cubic')

    cs = ax.contourf(x_gr,y_gr,z_gr, contour_levels, cmap=pl.cm.jet,
                             levels = np.linspace(cut_min,cut_max,contour_levels), extend='both') # , norm = LogNorm())
    #cs = ax.pcolor(x_gr,y_gr,z_gr, cmap=pl.cm.jet, vmin=cut_min, vmax=cut_max)
    
    cs.cmap.set_under('None')
    cs.cmap.set_over('k')

    ## dim1
    bx = 0.0
    by = lx[1] + ysep
    xls = np.linspace(bx,bx+lx[0],Nx[0])
    yls = np.linspace(by,by+lx[2],Nx[2])
    x_gr, y_gr = np.meshgrid(xls, yls)
    
    # grid the data.
    z_gr = griddata((x[1]+bx,y[1]+by), Zscore[1], (x_gr, y_gr), method='cubic')

    cs = ax.contourf(x_gr,y_gr,z_gr, contour_levels, cmap=pl.cm.jet,
                             levels = np.linspace(cut_min,cut_max,contour_levels), extend='both') # , norm = LogNorm())
    #cs = ax.pcolor(x_gr,y_gr,z_gr, cmap=pl.cm.jet, vmin=cut_min, vmax=cut_max)
    
    cs.cmap.set_under('None')
    cs.cmap.set_over('k')
    
    ## dim2
    bx = 0.0
    by = 0.0
    xls = np.linspace(bx,bx+lx[0],Nx[0])
    yls = np.linspace(by,by+lx[1],Nx[1])
    x_gr, y_gr = np.meshgrid(xls, yls)
    
    # grid the data.
    z_gr = griddata((x[2]+bx,y[2]+by), Zscore[2], (x_gr, y_gr), method='cubic')

    cs = ax.contourf(x_gr,y_gr,z_gr, contour_levels, cmap=pl.cm.jet,
                             levels = np.linspace(cut_min,cut_max,contour_levels), extend='both') # , norm = LogNorm())
    #cs = ax.pcolor(x_gr,y_gr,z_gr, cmap=pl.cm.jet, vmin=cut_min, vmax=cut_max)
    
    cs.cmap.set_under('None')
    cs.cmap.set_over('k')
    
    ## section  lines
    ax.plot([xsep,lx[0]],[sec[1]*dx[1], sec[1]*dx[1]],'k-',lw=0.5)
    ax.plot([xsep,lx[0] + xsep + lx[1]],[lx[1] + ysep + sec[2]*dx[2], lx[1] + ysep + sec[2]*dx[2]],'k-',lw=0.5)
    ax.plot([sec[0]*dx[0],sec[0]*dx[0]],[ysep,lx[1] + lx[2] + ysep],'k-',lw=0.5)
    ax.plot([lx[0] + xsep + sec[1]*dx[1],lx[0] + xsep + sec[1]*dx[1]],[lx[1] + ysep,lx[1] + ysep + lx[2]],'k-',lw=0.5)

    # scale bar
    ax.plot([12 + 200,108 + 200],[lx[1] - 14.0*ysep,lx[1] - 14.0*ysep],'k-',lw=10.0)
    fig.text(0.67,0.495,r'$\rm 100 ~mm$')

    # ax.plot(x, y, 'ko', markersize=4)

    ax.set_xlim(0,lx[0] + lx[1] + xsep)
    ax.set_ylim(0,lx[1] + lx[2] + ysep)

    #ax.tick_params(axis='x',which='minor',bottom='on')
    #ax.tick_params(axis='y',which='minor',bottom='on')

    ax.xaxis.set_minor_locator(mlx)
    ax.yaxis.set_minor_locator(mly)

    tit = r'$\rm time = %s ~(day)$' % realtime
    ax.set_xticks([])
    ax.set_yticks([])

    #if (i_y == 1):
    ax.set_title(tit)
    #ax.set_xticklabels([])
    #if (i_y == 0):
    #ax.set_xlabel(r'$\left( \phi - \pi \right) \sin{\theta} + \pi$')
    #ax.set_xticklabels([r'$\rm \pi/2$',r'$\rm \pi$',r'$\rm 3\pi/2$',r'$\rm 2 \pi$'])

    #ax.set_xlabel(labels[0])
    #ax.set_ylabel(labels[1])
    
    #ax.set_yticklabels([r'$\rm \pi/2$',r'$\rm \pi$'])
    #else:
    #pl.setp(ax.get_yticklabels(), visible=False)
    #ax.set_yticklabels([])

    #if (i_x == n_X - 1 and i_y == 0) :
    cbar_ax = fig.add_axes([0.93,0.15,0.02,0.7])
    cbar = pl.colorbar(cs,cax=cbar_ax)
    cbar.set_ticks([cut_min,cut_max])
    #mticks = cbar.norm([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20])
    #cbar.ax.yaxis.set_ticks(mticks, minor=True)
    #cbar.set_ticklabels([0, r'%.2e' % np.max(z)])
    cbar.ax.set_yticklabels([r'%.2e' % cut_min, r'%.2e' % cut_max],rotation=90)
    fig.text(0.96,0.5,lbl, rotation='vertical')
    
    ax.set_aspect('equal')

    # inset
    ax2 = fig.add_axes([0.52, 0.13, 0.38, 0.28], facecolor=(1.,1.,1.))
    ax2.hist(z, histtype='stepfilled',bins=50,density=True, fc='#CCCCCC', alpha=0.5, edgecolor='black', linewidth=1.2)  # bins='auto'
    ax2.set_xlabel(ag_names[me])
    #ax2.set_ylabel(r'$\rm probability$')
    
    #ax2.set_xlim(cut_min,cut_max)
    #ax2.set_xticklabels([r'%.2e' % cut_min, r'%.2e' % cut_max])
    
    ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    pl.savefig(fw, format = 'png', dpi=100, orientation='landscape')
    pl.close()

if __name__ == '__main__':
    main()
