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

rc('font',**{'family':'sans-serif','serif':['Times'], 'weight': 'bold', 'size' : 36})
rc('text', usetex=True)
rcParams.update({'figure.autolayout': True})

command = "mkdir corr"
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

AGENTS = ("mic","neu","sAb","fAb","ast","phr","tau","typ","group", "atrophy")
#ag_tit = (r'$\rm Microglia$', r'$\rm Neurons$', r'$\rm [sA\beta]$',r'$\rm [fA\beta]$',r'$\rm Astrogliosis$',r'$\rm [\tau]$',r'$\rm Tissue$',r'$\rm VOI$',r'$\rm Atrophy$')
ag_tit = (r'$M$', r'$N$', r'$S$',r'$F$',r'$A$',r'$P$',r'$T$',r'$\rm Tissue$',r'$\rm VOI$',r'$D$')
ag_names = (r'$M ~{\rm (mL^{-1})}$', r'$N ~{\rm (mL^{-1})}$', r'$S ~{\rm (\mu M)}$',r'$F ~{\rm (\mu M)}$',r'$\rm Astrogliosis$',r'$P ~{\rm (\mu M)}$',r'$T ~{\rm (\mu M)}$', r'$\rm Tissue$',r'$\rm VOI$',r'$\rm Atrophy$')
ag_sgn = ('M','N','S','F','A','P','T','TYPE','VOI','atrophy')

me_mic = 0
me_neu = 1
me_sAb = 2
me_fAb = 3
me_ast = 4
me_phr = 5
me_tau = 6
me_typ = 7
me_group = 8
me_atrophy = 9

x = []
y = []
Npoints = []

agents = []
agents0 = []

Pearson = []
minimum = []
maximum = []

ag = []
tis = []

Nx = np.zeros(4, dtype=np.int32)
dx = np.zeros(3, dtype=np.float32)
lx = np.zeros(3, dtype=np.float32)

###################################################################
###################################################################
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

###################################################################
###################################################################
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
    
    # dim0
    data_sec.append(data[sec[0]-1,:,:,0,:])
    data_sec0.append(data0[sec[0]-1,:,:,0,:])
    Npoints.append(Nx[1] * Nx[2])
    x.append(np.zeros(Npoints[0], dtype=np.float32))
    y.append(np.zeros(Npoints[0], dtype=np.float32))
    agents.append([])
    agents0.append([])
    
    # dim1
    data_sec.append(data[:,sec[1]-1,:,0,:])
    data_sec0.append(data0[:,sec[1]-1,:,0,:])
    Npoints.append(Nx[0] * Nx[2])
    x.append(np.zeros(Npoints[1], dtype=np.float32))
    y.append(np.zeros(Npoints[1], dtype=np.float32))
    agents.append([])
    agents0.append([])
    
    # dim2
    data_sec.append(data[:,:,sec[2]-1,0,:])
    data_sec0.append(data0[:,:,sec[2]-1,0,:])
    Npoints.append(Nx[0] * Nx[1])
    x.append(np.zeros(Npoints[2], dtype=np.float32))
    y.append(np.zeros(Npoints[2], dtype=np.float32))
    agents.append([])
    agents0.append([])
    
    for i in range(3):
        for j in range(Nx[3]):
            agents[i].append(np.zeros(Npoints[i], dtype=np.float32))
            agents0[i].append(np.zeros(Npoints[i], dtype=np.float32))
    
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
                agents0[0][ag_id][c] = data_sec0[0][i,j,ag_id]
            
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
                agents0[1][ag_id][c] = data_sec0[1][i,j,ag_id]
            
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
                agents0[2][ag_id][c] = data_sec0[2][i,j,ag_id]
            
            c += 1
    
    #prune_data()
    
    Zscore = []
    for me in range(Nx[3]-2):
        Zscore.append([])
        z = data[:,:,:,0,me].flatten()
        SMALL = 1.e-6
        z[z <= SMALL] = 0
        Z = np.ma.masked_equal(z,0)
        z = Z.compressed() # get normal array with masked values removed
        
        mu = np.mean(z)
        sig = np.std(z)
        
        zs = (z-mu)/sig
        if (math.isnan(sig)):
            continue
        
        cut_min = np.min(zs)
        cut_max = np.max(zs)
        
        for i in range(3):
            Zscore[me].append( (agents[i][me] - mu)/sig )

        ag.append( (data[:,:,:,0,me].flatten() - mu) / sig )
    tis.append(data[:,:,:,0,-1].flatten())

    for me in range(Nx[3]-2):
        Pearson.append([])
        minimum.append([])
        maximum.append([])
        for me2 in range(0,Nx[3]-2):
            Pearson[me].append([])
            
            for i in range(3):
                Pearson[me][me2].append(Zscore[me][i] * Zscore[me2][i] / (len(Zscore[me][i])))
                Pearson[me][me2][i][agents[i][Nx[3]-2] <= 0] = np.nan
                
                X = Pearson[me][me2][i][~np.isnan(Pearson[me][me2][i])]
                
                if (i == 0):
                    cut_max = np.max(X)
                    cut_min = np.min(X)
                else:
                    if (cut_max < np.max(X)):
                        cut_max = np.max(X)
                    if (cut_min > np.min(X)):
                        cut_min = np.min(X)
            
            if (cut_min < 0.0):
                cut_min *= -1.0
            if (cut_max < 0.0):
                cut_max *= -1.0
            if (cut_min < cut_max):
                cut_max = cut_min
            else:
                cut_min = cut_max
            cut_min *= -1.0
            
            minimum[me].append(cut_min)
            maximum[me].append(cut_max)


    # dummy image generation - preset
    ########################
    me = 0
    me2 =1
    fw = "corr/corr_" + AGENTS[me] + "_" + AGENTS[me2] + "_s" + str(sec[0]) + "_" + str(sec[1]) + "_" + str(sec[2]) + "_t" + time
    #make_plot(me,me2,fw,100)
    ########################

    for me in range(Nx[3]-2):
        for me2 in range(me+1,Nx[3]-2):
            print("ploting", ag_sgn[me],ag_sgn[me2])
            fw = "corr/corr_" + AGENTS[me] + "_" + AGENTS[me2] + "_s" + str(sec[0]) + "_" + str(sec[1]) + "_" + str(sec[2]) + "_t" + time
            #make_plot(me,me2,fw,100)
    
    # combine all diagrams
    combine_plots()

###################################################################
###################################################################
def combine_plots():
    me_list = ( me_mic, me_neu, me_sAb, me_fAb, me_ast, me_tau )

    fig = pl.figure(figsize=(20,20))

    n_X = 6
    n_Y = 6
    
    b_X = 0.03
    b_Y = 0.03
    
    sep_X = 0.01
    sep_Y = 0.01
    
    l_X = (1.0 - sep_X * n_X - b_X ) / n_X
    l_Y = (1.0 - sep_Y * n_Y - b_Y ) / n_Y
    
    # from top left to bottom right
    for j in range(n_Y):
        y_ax = 1.0 - b_Y - l_Y - (l_Y + sep_Y) * j
        for i in range(n_X):
            x_ax = b_X + (l_X + sep_X) * i
            print("combine plot", ag_sgn[j],ag_sgn[i])
            
            if (i < j) :
                ax = fig.add_axes([x_ax, y_ax, l_X, l_Y], facecolor=(1.,1.,1.))
                cs = add_contour(fig,ax,j,i,100)
                cbar_ax = fig.add_axes([x_ax + 0.50*l_X, y_ax + 0.05 * l_Y, 0.04*l_X, 0.4*l_Y])
                cbar = pl.colorbar(cs,cax=cbar_ax)
                cbar.set_ticks([minimum[i][j],0,maximum[i][j]])
                cbar.ax.set_yticklabels([r'$\rm %.1e$' % minimum[i][j], '', r'$\rm %.1e$' % maximum[i][j]],rotation=0,fontsize=24)
                
            elif (i > j) :
                ax = fig.add_axes([x_ax + 0.05*l_X, y_ax + 0.05*l_Y, 0.9*l_X, 0.9*l_Y], facecolor=(1.,1.,1.))
                add_scatter(fig,ax,i,j)
            else :
                fig.text(x_ax + 0.5*l_X, y_ax + 0.5*l_Y, ag_tit[me_list[j]], ha='center', va='center', fontsize=40)
                #ax.set_axis_off()
                continue
    
    tit = r'$\rm t = %g ~(day)$' % realtime
    ax.set_xticks([])
    ax.set_yticks([])

    lbl = r'${\rm Scatter ~plots} ~z_i ~{\rm vs.} ~z_j$'
    fig.text(0.5,0.98,lbl, ha='center', va='center', fontsize=40)
    
    lbl = r'${\rm Correlation ~maps} ~{\rm corr}{(z_i, z_j)}$'
    fig.text(0.02,0.5,lbl, ha='center', va='center', rotation=90, fontsize=40)

    fw = "corr/corr_all_s" + str(sec[0]) + "_" + str(sec[1]) + "_" + str(sec[2]) + "_t" + time
    pl.savefig(fw + '.png', format = 'png', dpi=200, orientation='landscape')
    pl.close()

###################################################################
###################################################################
def add_contour(fig,ax,me,me2,contour_levels):
    P = Pearson[me][me2]
    cut_min = minimum[me][me2]
    cut_max = maximum[me][me2]

    xsep = 1.
    ysep = 1.
    
    mlx   = MultipleLocator(50)
    mly   = MultipleLocator(50)

    ## dim0
    bx = lx[0] + xsep
    by = lx[1] + ysep
    xls = np.linspace(bx,bx+lx[1],Nx[1])
    yls = np.linspace(by,by+lx[2],Nx[2])
    x_gr, y_gr = np.meshgrid(xls, yls)
    
    # grid the data.
    z_gr = griddata((x[0]+bx,y[0]+by), P[0], (x_gr, y_gr), method='nearest')
    
    cs = ax.contourf(x_gr,y_gr,z_gr, contour_levels, cmap=pl.cm.coolwarm,
                             levels = np.linspace(cut_min,cut_max,contour_levels), extend='both')

    ## dim1
    bx = 0.0
    by = lx[1] + ysep
    xls = np.linspace(bx,bx+lx[0],Nx[0])
    yls = np.linspace(by,by+lx[2],Nx[2])
    x_gr, y_gr = np.meshgrid(xls, yls)
    
    # grid the data.
    z_gr = griddata((x[1]+bx,y[1]+by), P[1], (x_gr, y_gr), method='nearest')

    cs = ax.contourf(x_gr,y_gr,z_gr, contour_levels, cmap=pl.cm.coolwarm,
                             levels = np.linspace(cut_min,cut_max,contour_levels), extend='both')
    
    ## dim2
    bx = 0.0
    by = 0.0
    xls = np.linspace(bx,bx+lx[0],Nx[0])
    yls = np.linspace(by,by+lx[1],Nx[1])
    x_gr, y_gr = np.meshgrid(xls, yls)
    
    # grid the data.
    z_gr = griddata((x[2]+bx,y[2]+by), P[2], (x_gr, y_gr), method='nearest')

    cs = ax.contourf(x_gr,y_gr,z_gr, contour_levels, cmap=pl.cm.coolwarm,
                             levels = np.linspace(cut_min,cut_max,contour_levels), extend='both')
    
    ## section  lines
    ax.plot([xsep,lx[0]],[sec[1]*dx[1], sec[1]*dx[1]],'k-',lw=0.5)
    ax.plot([xsep,lx[0] + xsep + lx[1]],[lx[1] + ysep + sec[2]*dx[2], lx[1] + ysep + sec[2]*dx[2]],'k-',lw=0.5)
    ax.plot([sec[0]*dx[0],sec[0]*dx[0]],[ysep,lx[1] + lx[2] + ysep],'k-',lw=0.5)
    ax.plot([lx[0] + xsep + sec[1]*dx[1],lx[0] + xsep + sec[1]*dx[1]],[lx[1] + ysep,lx[1] + ysep + lx[2]],'k-',lw=0.5)

    # scale bar
    if (me == 0):
        from matplotlib.patches import Rectangle
        ax.add_patch(Rectangle((150, 205), 100, 2, alpha=1, color='k'))
        ax.text(200,150,r'$\rm 100 ~mm$')

    ax.set_xlim(0,lx[0] + lx[1] + xsep)
    ax.set_ylim(0,lx[1] + lx[2] + ysep)

    ax.xaxis.set_minor_locator(mlx)
    ax.yaxis.set_minor_locator(mly)
    
    ax.set_xticks([])
    ax.set_yticks([])

    ax.set_aspect('equal')
    
    return cs

###################################################################
###################################################################
def add_scatter(fig,ax,me,me2):
    marray = np.ma.masked_where(tis[0] <= 0,tis[0]).mask
    z11 = np.ma.array(ag[me], mask = marray).compressed()
    z12 = np.ma.array(ag[me2], mask = marray).compressed()
    ax.hist2d(z11, z12, bins=200, cmap='Blues',alpha=1, norm = matplotlib.colors.LogNorm())
    
    if (me - 1 == me2) :
        ax.tick_params(axis='both', which='major', labelsize=20)
        lblx = r'$z_{%s}$' % ag_sgn[me]
        ax.set_xlabel(lblx, fontsize = 32)
        lbly = r'$z_{%s}$' % ag_sgn[me2]
        ax.set_ylabel(lbly, fontsize = 32)
    else :
        ax.set_xticks([])
        ax.set_yticks([])
    
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(1))

###################################################################
###################################################################
def prune_data():
    # find useless data
    for i in range(3):
        marray = np.ma.masked_where(agents[i][-1] < 0,agents[i][-1]).mask
        
        x[i] = np.ma.array(x[i], mask = marray).compressed()
        y[i] = np.ma.array(y[i], mask = marray).compressed()
    
        for j in range(Nx[3]):
            agents[i][j] = np.ma.array(agents[i][j], mask = marray).compressed()

###################################################################
###################################################################
def make_plot(me,me2,fw,contour_levels):
    
    P = Pearson[me][me2]
    cut_min = minimum[me][me2]
    cut_max = maximum[me][me2]
    
    fig = pl.figure(figsize=(9,9))
    #ax = fig.add_subplot(111)
    ax = fig.add_axes([0.02,0.02,0.9,0.9])

    xsep = 1.
    ysep = 1.
    
    mlx   = MultipleLocator(10)
    mly   = MultipleLocator(10)

    lbl = r'${\rm corr}(%s,%s)$' % (ag_sgn[me],ag_sgn[me2])

    ## dim0
    bx = lx[0] + xsep
    by = lx[1] + ysep
    xls = np.linspace(bx,bx+lx[1],Nx[1])
    yls = np.linspace(by,by+lx[2],Nx[2])
    x_gr, y_gr = np.meshgrid(xls, yls)
    
    # grid the data.
    z_gr = griddata((x[0]+bx,y[0]+by), P[0], (x_gr, y_gr), method='nearest')
    
    cs = ax.contourf(x_gr,y_gr,z_gr, contour_levels, cmap=pl.cm.coolwarm,
                             levels = np.linspace(cut_min,cut_max,contour_levels), extend='both')
    #cs = ax.pcolor(x_gr,y_gr,z_gr, cmap=pl.cm.coolwarm, vmin=cut_min, vmax=cut_max)
    
    #cs.cmap.set_under('None')
    #cs.cmap.set_over('None')

    ## dim1
    bx = 0.0
    by = lx[1] + ysep
    xls = np.linspace(bx,bx+lx[0],Nx[0])
    yls = np.linspace(by,by+lx[2],Nx[2])
    x_gr, y_gr = np.meshgrid(xls, yls)
    
    # grid the data.
    z_gr = griddata((x[1]+bx,y[1]+by), P[1], (x_gr, y_gr), method='nearest')

    cs = ax.contourf(x_gr,y_gr,z_gr, contour_levels, cmap=pl.cm.coolwarm,
                             levels = np.linspace(cut_min,cut_max,contour_levels), extend='both') # , norm = LogNorm())
    #cs = ax.pcolor(x_gr,y_gr,z_gr, cmap=pl.cm.coolwarm, vmin=cut_min, vmax=cut_max)
    
    #cs.cmap.set_under('None')
    #cs.cmap.set_over('None')
    
    ## dim2
    bx = 0.0
    by = 0.0
    xls = np.linspace(bx,bx+lx[0],Nx[0])
    yls = np.linspace(by,by+lx[1],Nx[1])
    x_gr, y_gr = np.meshgrid(xls, yls)
    
    # grid the data.
    z_gr = griddata((x[2]+bx,y[2]+by), P[2], (x_gr, y_gr), method='nearest')

    cs = ax.contourf(x_gr,y_gr,z_gr, contour_levels, cmap=pl.cm.coolwarm,
                             levels = np.linspace(cut_min,cut_max,contour_levels), extend='both') # , norm = LogNorm())
    #cs = ax.pcolor(x_gr,y_gr,z_gr, cmap=pl.cm.coolwarm, vmin=cut_min, vmax=cut_max)
    
    #cs.cmap.set_under('None')
    #cs.cmap.set_over('None')
    
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

    tit = r'$\rm corr($%s$, $%s$) \rm , ~t = %g ~(day)$' % (ag_tit[me],ag_tit[me2],realtime)
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
    cbar_ax = fig.add_axes([0.93,0.20,0.02,0.7])
    cbar = pl.colorbar(cs,cax=cbar_ax)
    cbar.set_ticks([cut_min,0,cut_max])
    #mticks = cbar.norm([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20])
    #cbar.ax.yaxis.set_ticks(mticks, minor=True)
    #cbar.set_ticklabels([0, r'%.2e' % np.max(z)])
    cbar.ax.set_yticklabels([r'$\rm %.1e$' % cut_min, '', r'$\rm %.1e$' % cut_max],rotation=90,fontsize=30)
    #cbar.ax.set_yticklabels([r'$\rm %.1f$' % cut_min, r'$\rm %.1f$' % cut_max],rotation=90,fontsize=30)
    fig.text(0.96,0.5,lbl, rotation='vertical', fontsize=30)
    
    ax.set_aspect('equal')

    # inset
    ax2 = fig.add_axes([0.56, 0.11, 0.32, 0.32], facecolor=(1.,1.,1.))
    marray = np.ma.masked_where(tis[0] <= 0,tis[0]).mask
    z11 = np.ma.array(ag[me], mask = marray).compressed()
    z12 = np.ma.array(ag[me2], mask = marray).compressed()
    ax2.hist2d(z11, z12, bins=200, cmap='Blues',alpha=1, norm = matplotlib.colors.LogNorm())
    
    ax2.tick_params(axis='both', which='major', labelsize=20)

    lblx = r'$z_{%s}$' % ag_sgn[me]
    ax2.set_xlabel(lblx, fontsize = 26)
    lbly = r'$z_{%s}$' % ag_sgn[me2]
    ax2.set_ylabel(lbly, fontsize = 26)

    pl.savefig(fw + '.png', format = 'png', dpi=100, orientation='landscape')
    pl.close()

###################################################################
###################################################################
if __name__ == '__main__':
    main()
