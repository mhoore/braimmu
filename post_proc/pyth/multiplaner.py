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
import matplotlib.image as mpimg

rc('font',**{'family':'sans-serif','serif':['Times'], 'weight': 'bold', 'size' : 26})
rc('text', usetex=True)
rcParams.update({'figure.autolayout': True})

command = "mkdir multi"
p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

# Data Input
scale_fac = 1.e-3

if len(sys.argv) < 8:
    print(sys.argv)
    sys.exit("Error: chack the arguments")

sec = []
fr = sys.argv[1]
sec.append(int(sys.argv[2]))
sec.append(int(sys.argv[3]))
sec.append(int(sys.argv[4]))
time = sys.argv[5]
realtime = float(sys.argv[6])
fr0 = sys.argv[7]

AGENTS = ("mic","neu","sAb","fAb","ast","typ","group", "atrophy")
ag_tit = (r'$\rm Microglia$', r'$\rm Neurons$', r'$\rm [sA\beta]$',r'$\rm [fA\beta]$',r'$\rm Astrogliosis$',r'$\rm Tissue$',r'$\rm VOI$',r'$\rm Atrophy$')
ag_names = (r'$M ~{\rm (mL^{-1})}$', r'$N ~{\rm (mL^{-1})}$', r'$S ~{\rm (\mu M)}$',r'$F ~{\rm (\mu M)}$',r'$\rm Astrogliosis$',r'$\rm Tissue$',r'$\rm VOI$',r'$\rm Atrophy$')
ag_sgn = ('M','N','S','F','A','TYPE','VOI','atrophy')

me_mic = 0
me_neu = 1
me_sAb = 2
me_fAb = 3
me_ast = 4
me_typ = 5
me_group = 6
me_atrophy = 7

class Volume_of_Interest:
    # Constructor
    def __init__(self, label, name, descr):
        self.label = label
        self.name = name
        self.descr = descr

### setup the VOIs
data = np.genfromtxt("/home/mho18/braimmu/base/voi.txt", delimiter=",", dtype="|U40", comments='#')

voi = []
for i in range(len(data)):
    voi.append( Volume_of_Interest(int(data[i][0]) , data[i][1], data[i][2]) )
data = []

random_color = np.random.rand (len(voi),3)

x = []
y = []
Npoints = []

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

    img0 = nib.load(fr0)
    data0 = img0.get_fdata()
    
    Nx[0], Nx[1], Nx[2], dum, Nx[3] = data.shape
    dx[0], dx[1], dx[2], dt, dum = img.header.get_zooms()
    
    for i in range(3):
        dx[i] *= scale_fac
        lx[i] = Nx[i] * dx[i]
    
    # dim0
    Npoints.append(Nx[1] * Nx[2])
    x.append(np.zeros(Npoints[0], dtype=np.float32))
    y.append(np.zeros(Npoints[0], dtype=np.float32))
    
    # dim1
    Npoints.append(Nx[0] * Nx[2])
    x.append(np.zeros(Npoints[1], dtype=np.float32))
    y.append(np.zeros(Npoints[1], dtype=np.float32))
    
    # dim2
    Npoints.append(Nx[0] * Nx[1])
    x.append(np.zeros(Npoints[2], dtype=np.float32))
    y.append(np.zeros(Npoints[2], dtype=np.float32))
    
    # dim 0
    ii = 1
    jj = 2
    c = 0
    for i in range(Nx[ii]):
        dum = dx[ii] * i
        for j in range(Nx[jj]):
            x[0][c] = dum
            y[0][c] = dx[jj] * j
            
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
            
            c += 1
    
    #prune_data()
    
    # dummy image generation - preset
    ########################
    me = 0
    z = data[:,:,:,0,me].flatten()
    SMALL = 1.e-6
    z[z <= SMALL] = 0
    Z = np.ma.masked_equal(z,0)
    z = Z.compressed() # get normal array with masked values removed

    mu = np.mean(z)
    sig = np.std(z)

    zs = (z-mu)/sig
    if (math.isnan(sig) == False):
        cut_min = np.min(zs)
        cut_max = np.max(zs)

        Zscore = []
        Zscore.append( ( data[sec[0]-1,:,:,0,me].flatten() - mu)/sig )
        Zscore.append( ( data[:,sec[1]-1,:,0,me].flatten() - mu)/sig )
        Zscore.append( ( data[:,:,sec[2]-1,0,me].flatten() - mu)/sig )

        fw = "multi/m_Z_" + AGENTS[me] + "_s" + str(sec[0]) + "_" + str(sec[1]) + "_" + str(sec[2]) + "_t" + time

        z   = data[:,:,:,0,me].flatten()
        tis = data[:,:,:,0,-2].flatten()

        make_plot(Zscore,z,tis,me,cut_min,cut_max,fw,100)
    ########################

    for me in range(Nx[3]-2):
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
        
        Zscore = []
        Zscore.append( ( data[sec[0]-1,:,:,0,me].flatten() - mu)/sig )
        Zscore.append( ( data[:,sec[1]-1,:,0,me].flatten() - mu)/sig )
        Zscore.append( ( data[:,:,sec[2]-1,0,me].flatten() - mu)/sig )
        
        fw = "multi/m_Z_" + AGENTS[me] + "_s" + str(sec[0]) + "_" + str(sec[1]) + "_" + str(sec[2]) + "_t" + time

        z   = data[:,:,:,0,me].flatten()
        tis = data[:,:,:,0,-2].flatten()

        make_plot(Zscore,z,tis,me,cut_min,cut_max,fw,100)

    # tissue
    me = Nx[3]-2
    z = data[:,:,:,0,me].flatten()
    Zscore = []
    Zscore.append( data[sec[0]-1,:,:,0,me].flatten() )
    Zscore.append( data[:,sec[1]-1,:,0,me].flatten() )
    Zscore.append( data[:,:,sec[2]-1,0,me].flatten() )
    
    fw = "multi/m_type_s" + str(sec[0]) + "_" + str(sec[1]) + "_" + str(sec[2]) + "_t" + time
    make_plot(Zscore,z,tis,me,-1,3,fw,4)

    # voi
    me = Nx[3]-1
    
    data_raw = data[:,:,:,0,me]
    data_voi = np.zeros_like(data_raw)
    for i in range(len(voi)):
        marray = np.ma.masked_where(data_raw - voi[i].label == 0, data_raw)
        data_voi += marray.mask * (i + 1)
    
    z = data_voi.flatten()
    Zscore = []
    Zscore.append( data_voi[sec[0]-1,:,:].flatten() )
    Zscore.append( data_voi[:,sec[1]-1,:].flatten() )
    Zscore.append( data_voi[:,:,sec[2]-1].flatten() )
    fw = "multi/m_voi_s" + str(sec[0]) + "_" + str(sec[1]) + "_" + str(sec[2]) + "_t" + time
    make_plot(Zscore,z,tis,me,0,len(voi),fw,len(voi))
    
    # attrophy (1-N)
    me = 1
    atrophy = []
    atrophy.append(- data[sec[0]-1,:,:,0,me].flatten() / data0[sec[0]-1,:,:,0,me].flatten() + 1.0)
    atrophy.append(- data[:,sec[1]-1,:,0,me].flatten() / data0[:,sec[1]-1,:,0,me].flatten() + 1.0)
    atrophy.append(- data[:,:,sec[2]-1,0,me].flatten() / data0[:,:,sec[2]-1,0,me].flatten() + 1.0)
    
    z = - data[:,:,:,0,me].flatten() / data0[:,:,:,0,me].flatten() + 1.0
    SMALL = 1.e-6
    z[z <= SMALL] = 0
    Z = np.ma.masked_equal(z,0)
    z = Z.compressed() # get normal array with masked values removed
    
    cut_min = 0.0
    cut_max = 1.0
    
    fw = "multi/m_atrophy_s" + str(sec[0]) + "_" + str(sec[1]) + "_" + str(sec[2]) + "_t" + time

    z = - data[:,:,:,0,me].flatten() / data0[:,:,:,0,me].flatten() + 1.0
    tis = data[:,:,:,0,-2].flatten()

    make_plot(atrophy,z,tis,Nx[3],cut_min,cut_max,fw,100)
    
    data = []
    data_sec = []
    data0 = []
    data_sec0 = []
    
    print('HERE')
    # combine all diagrams
    combine_plots()

## plots
def combine_plots():
    me_1 = ( me_mic, me_neu, me_atrophy, me_fAb, me_sAb, me_ast )
    tit  = (   r'a',   r'b',       r'c',   r'd',   r'e',   r'f' )
    
    fig = pl.figure(figsize=(15,10))

    n_X = 3
    n_Y = 2
    
    b_X = 0.01
    b_Y = 0.01
    
    sep_X = 0.01
    sep_Y = 0.01
    
    l_X = 0.96 / n_X
    l_Y = 0.97 / n_Y
    
    c= 0
    # from top left to bottom right
    for j in range(n_Y):
        y_ax = 1.0 - b_Y - l_Y - (l_Y + sep_Y) * j
        for i in range(n_X):
            x_ax = b_X + (l_X + sep_X) * i
            
            ax = fig.add_axes([x_ax, y_ax, l_X, l_Y], facecolor=(1.,1.,1.))
            ax.set_axis_off()
            
            fw = "multi/m_Z_" + AGENTS[me_1[c]] + "_s" + str(sec[0]) + "_" + str(sec[1]) + "_" + str(sec[2]) + "_t" + time + ".png"
            if (me_1[c] == me_atrophy):
                fw = "multi/m_" + AGENTS[me_1[c]] + "_s" + str(sec[0]) + "_" + str(sec[1]) + "_" + str(sec[2]) + "_t" + time + ".png"
            im = mpimg.imread(fw)
            ax.imshow(im)
            
            fig.text(x_ax + 0.01, y_ax + l_Y - 0.03, tit[c], fontsize=30)
            
            c += 1

    fw = "multi/m_all_s" + str(sec[0]) + "_" + str(sec[1]) + "_" + str(sec[2]) + "_t" + time
    pl.savefig(fw + '.png', format = 'png', dpi=200, orientation='landscape')
    pl.close()
    
"""
def prune_data():
    # find useless data
    for i in range(3):
        marray = np.ma.masked_where(data[sec[0]-1,:,:,0,-1] < 0,data[sec[0]-1,:,:,0,-1]).mask
        
        x[i] = np.ma.array(x[i], mask = marray).compressed()
        y[i] = np.ma.array(y[i], mask = marray).compressed()
    
        for j in range(Nx[3]):
            agents[i][j] = np.ma.array(agents[i][j], mask = marray).compressed()
"""
## plots
def make_plot(Zscore,z,tis,me,cut_min,cut_max,fw,contour_levels):
    
    fig = pl.figure(figsize=(9,9))
    #ax = fig.add_subplot(111)
    ax = fig.add_axes([0.02,0.02,0.9,0.9])

    xsep = 1.
    ysep = 1.
    
    mlx   = MultipleLocator(10)
    mly   = MultipleLocator(10)

    cmap = pl.cm.jet
    lbl = r'$z_{%s}$' % ag_sgn[me]
    if (me == Nx[3]-2):
        lbl = r'$\rm Tissue$'
    elif (me == Nx[3]-1):
        lbl = r'$\rm VOI$'
        cmap = matplotlib.colors.ListedColormap (random_color)
        #cmap = pl.cm.get_cmap('Greys', len(voi))
    elif (me == Nx[3]):
        lbl = r'$\rm Atrophy$'


    ## dim0
    bx = lx[0] + xsep
    by = lx[1] + ysep
    xls = np.linspace(bx,bx+lx[1],Nx[1])
    yls = np.linspace(by,by+lx[2],Nx[2])
    x_gr, y_gr = np.meshgrid(xls, yls)
    
    # grid the data.
    z_gr = griddata((x[0]+bx,y[0]+by), Zscore[0], (x_gr, y_gr), method='nearest')

    print(np.min(z_gr),cut_min)
    cs = ax.contourf(x_gr,y_gr,z_gr, contour_levels, cmap=cmap,
                             levels = np.linspace(cut_min,cut_max,contour_levels), extend='both') # , norm = LogNorm())
    #cs = ax.pcolor(x_gr,y_gr,z_gr, cmap=cmap, vmin=cut_min, vmax=cut_max)
    
    cs.cmap.set_under('None')
    cs.cmap.set_over('k')

    ## dim1
    bx = 0.0
    by = lx[1] + ysep
    xls = np.linspace(bx,bx+lx[0],Nx[0])
    yls = np.linspace(by,by+lx[2],Nx[2])
    x_gr, y_gr = np.meshgrid(xls, yls)
    
    # grid the data.
    z_gr = griddata((x[1]+bx,y[1]+by), Zscore[1], (x_gr, y_gr), method='nearest')

    cs = ax.contourf(x_gr,y_gr,z_gr, contour_levels, cmap=cmap,
                             levels = np.linspace(cut_min,cut_max,contour_levels), extend='both') # , norm = LogNorm())
    #cs = ax.pcolor(x_gr,y_gr,z_gr, cmap=cmap, vmin=cut_min, vmax=cut_max)
    
    cs.cmap.set_under('None')
    cs.cmap.set_over('k')
    
    ## dim2
    bx = 0.0
    by = 0.0
    xls = np.linspace(bx,bx+lx[0],Nx[0])
    yls = np.linspace(by,by+lx[1],Nx[1])
    x_gr, y_gr = np.meshgrid(xls, yls)
    
    # grid the data.
    z_gr = griddata((x[2]+bx,y[2]+by), Zscore[2], (x_gr, y_gr), method='nearest')

    cs = ax.contourf(x_gr,y_gr,z_gr, contour_levels, cmap=cmap,
                             levels = np.linspace(cut_min,cut_max,contour_levels), extend='both') # , norm = LogNorm())
    #cs = ax.pcolor(x_gr,y_gr,z_gr, cmap=cmap, vmin=cut_min, vmax=cut_max)
    
    cs.cmap.set_under('None')
    cs.cmap.set_over('None')
    
    ## section  lines
    ax.plot([xsep,lx[0]],[sec[1]*dx[1], sec[1]*dx[1]],'k-',lw=0.5)
    ax.plot([xsep,lx[0] + xsep + lx[1]],[lx[1] + ysep + sec[2]*dx[2], lx[1] + ysep + sec[2]*dx[2]],'k-',lw=0.5)
    ax.plot([sec[0]*dx[0],sec[0]*dx[0]],[ysep,lx[1] + lx[2] + ysep],'k-',lw=0.5)
    ax.plot([lx[0] + xsep + sec[1]*dx[1],lx[0] + xsep + sec[1]*dx[1]],[lx[1] + ysep,lx[1] + ysep + lx[2]],'k-',lw=0.5)

    # scale bar
    from matplotlib.patches import Rectangle
    #cax = pl.gca()
    ax.add_patch(Rectangle((150, 205), 100, 2, alpha=1, color='k'))
    #ax.plot([12 + 200,108 + 200],[lx[1] - 14.0*ysep,lx[1] - 14.0*ysep],'k-',lw=10.0)
    fig.text(0.39,0.47,r'$\rm 100 ~mm$')

    # ax.plot(x, y, 'ko', markersize=4)

    ax.set_xlim(0,lx[0] + lx[1] + xsep)
    ax.set_ylim(0,lx[1] + lx[2] + ysep)

    #ax.tick_params(axis='x',which='minor',bottom='on')
    #ax.tick_params(axis='y',which='minor',bottom='on')

    ax.xaxis.set_minor_locator(mlx)
    ax.yaxis.set_minor_locator(mly)

    tit = r'%s ,  $\rm t = %g ~(day)$' % (ag_tit[me],realtime)
    ax.set_xticks([])
    ax.set_yticks([])

    #if (i_y == 1):
    fig.text(0.50,0.96,tit,fontsize=30, ha='center', va='center')
    #ax.set_title(tit)
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
    #cbar.ax.set_yticklabels([r'$\rm %.1e$' % cut_min, r'$\rm %.1e$' % cut_max],rotation=90,fontsize=30)
    cbar.ax.set_yticklabels([r'$\rm %.1f$' % cut_min, r'$\rm %.1f$' % cut_max],rotation=90,fontsize=30)
    fig.text(0.96,0.5,lbl, rotation='vertical', fontsize=30)
    
    ax.set_aspect('equal')

    # inset
    ax2 = fig.add_axes([0.48, 0.11, 0.38, 0.28], facecolor=(1.,1.,1.))
    if (me == Nx[3]-2): # for type (tissue) histogram
        Z = np.ma.masked_equal(z,-1)
        z = Z.compressed()
        ax2.hist(z, histtype='stepfilled',bins=50,density=False, fc='grey', alpha=0.6, edgecolor='black', linewidth=1.2)  # bins='auto'
        ax2.set_xticklabels([r'$\rm -1 (empty)$',r'$\rm 0 (CSF)$', r'$\rm 1 (WM)$', r'$\rm 2 (GM)$'])
    elif (me == Nx[3]-1): # for group (voi) histogram
        Z = np.ma.masked_equal(z,0)
        z = Z.compressed()
        ax2.hist(z, histtype='stepfilled',bins=256,density=False, fc='grey', alpha=0.6, edgecolor='black', linewidth=1.2)  # bins='auto'
    else:
        marray = np.ma.masked_where(tis != 1,tis).mask
        z1 = np.ma.array(z, mask = marray).compressed()
        marray = np.ma.masked_where(tis != 2,tis).mask
        z2 = np.ma.array(z, mask = marray).compressed()
        
        SMALL = 1.e-6
        
        z1[z1 <= SMALL] = 0
        Z = np.ma.masked_equal(z1,0)
        z1 = Z.compressed()

        z2[z2 <= SMALL] = 0
        Z = np.ma.masked_equal(z2,0)
        z2 = Z.compressed()
        
        ax2.hist(z2, histtype='stepfilled',bins=50,density=False, fc='grey', alpha=0.6, edgecolor='black', linewidth=1.2, label=r'$\rm GM$')  # bins='auto'
        ax2.hist(z1, histtype='stepfilled',bins=50,density=False, fc='yellow', alpha=0.8, edgecolor='black', linewidth=1.2, label=r'$\rm WM$')  # bins='auto'
        ax2.legend(loc='upper right', fontsize=18, ncol=2)

        ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        #ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        
        a = np.histogram(z1,bins=50,density=False)
        b = np.histogram(z2,bins=50,density=False)
        
        ax2.set_ylim(0,1.30*max(np.max(a[0]),np.max(b[0])))
        
    #ax2.set_xlim(cut_min,cut_max)
    ax2.set_yticklabels([])
    
    ax2.set_title(ag_names[me], fontsize = 26)
    #ax2.set_xlabel(ag_names[me], fontsize = 26)
    #ax2.set_ylabel(r'$\rm probability$')

    pl.savefig(fw + '.png', format = 'png', dpi=100, orientation='landscape')
    pl.close()

if __name__ == '__main__':
    main()
