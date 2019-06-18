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

rc('font',**{'family':'sans-serif','serif':['Times'], 'weight': 'bold', 'size' : 20})
rc('text', usetex=True)
rcParams.update({'figure.autolayout': True})

# Data Input
scale_fac = 1.e-3
scale_time = 1.e-3 / 365.

if len(sys.argv) < 2:
    print(sys.argv)
    sys.exit("Error: chack the arguments")

fname = sys.argv[1]

AGENTS = ("mic","neu","sAb","fAb","ast","typ","group")
ag_tit = (r'$M/M_{\rm o}$', r'$N/N_{\rm max}$', r'$S / S_{\rm max}$',r'$F / F_{\rm max}$',r'$A$',r'$\rm Tissue$',r'$\rm VOI$',r'$\rm Atrophy$')
ag_names = (r'$M ~{\rm (mL^{-1})}$', r'$N ~{\rm (mL^{-1})}$', r'$S ~{\rm (\mu M)}$',r'$F ~{\rm (\mu M)}$',r'$\rm Astrogliosis$',r'$\rm Tissue$',r'$\rm VOI$',r'$\rm Atrophy$')
ag_sgn = ('M','N','S','F','A','TYPE','VOI','atrophy')

VOIs = []
VOIs.append( (r'$\rm FrW$', 30, 17,    'FLW ', 'FRW ', 'black',         2.2, '-' ) )
VOIs.append( (r'$\rm FrG$', 210, 211,  'FLG ', 'FRG ', 'blue',          2.2, '--' ) )
VOIs.append( (r'$\rm TeW$', 83, 59,    'TLW ', 'TRW ', 'red',           2.0, '-' ) )
VOIs.append( (r'$\rm TeG$', 218, 219,  'TLG ', 'TRG ', 'green',         2.0, '--' ) )
VOIs.append( (r'$\rm PaW$', 57, 105,   'PLW ', 'PRW ', 'cyan',          1.9, '-' ) )
VOIs.append( (r'$\rm PaG$', 6, 2,      'PLG ', 'PRG ', 'magenta',       1.9, '--' ) )
VOIs.append( (r'$\rm OcW$', 73, 45,    'OLW ', 'ORW ', 'olive',         1.8, '-' ) )
VOIs.append( (r'$\rm OcG$', 8, 4,      'OLG ', 'ORG ', 'lightgreen',    1.8, '--' ) )
VOIs.append( (r'$\rm Cer$', 67, 76,    'CeL ', 'CeR ', 'firebrick',     1.7, '-' ) )
VOIs.append( (r'$\rm Cau$', 39, 53,    'CaL ', 'CaR ', 'orange',        1.7, '--' ) )
VOIs.append( (r'$\rm Put$', 14, 16,    'PuL ', 'PuR ', 'steelblue',     1.6, '-' ) )
VOIs.append( (r'$\rm Tha$', 102, 203,  'ThL ', 'ThR ', 'hotpink',       1.6, '--' ) )
VOIs.append( (r'$\rm sTN$', 33, 23,    'sTNL', 'sTNR', 'slategray',     1.5, '-' ) )
VOIs.append( (r'$\rm GlP$', 12, 11,    'GPL ', 'GPR ', 'blueviolet',    1.5, '--' ) )

#VOIs.append( (r'$\rm LaV $', 3, 9,      'LVL ', 'LVR ', 'green',   2.2, '-' ) )
#VOIs.append( (r'$\rm Frx $', 29, 254,   'FoL ', 'FoR ', 'green',   1.8, '-' ) )

voi_id = np.zeros(len(VOIs),dtype=int)
voi_tit = np.zeros(len(VOIs),dtype='|U10')
voi_col = np.zeros(len(VOIs),dtype='|U10')

class Subplots:
    """ Arrange subplot grid structure (square box is assumed)"""

    totcnt = -1             # Total number of subplots

    # Constructor
    def __init__(self, f, lx, ly, sx, sy, bx, by, t):
        self.fig = f                    # Figure axes handle
        self.length_x = lx      # Length of the subplot box
        self.length_y = ly      # Length of the subplot box
        self.sepx = sx                  # Separation distance between subplots
        self.sepy = sy                  # Separation distance between subplots

        self.begx = bx                  # Beginning (offset) in the figure box
        self.begy = by                  # Beginning (offset) in the figure box
        self.tot = t                    # Total number of subplots in x direction

    # Add a subplot in the grid structure
    def addSubplot(self):

        # Increase the number of subplots in the figure
        self.totcnt += 1

        # Indices of the subplot in the figure
        self.nx = self.totcnt % (self.tot)
        self.ny = int(self.totcnt / (self.tot))

        self.xbeg = self.begx + self.nx*self.length_x + self.nx*self.sepx
        self.ybeg = self.begy + self.ny*self.length_y + self.ny*self.sepy
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length_x,self.length_y])

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

x = []
y = []

Nx = np.zeros(4, dtype=np.int32)
dx = np.zeros(3, dtype=np.float32)
lx = np.zeros(3, dtype=np.float32)

def extract(data_mn, data_sd, vlabel, VOIsL, VOIsR):
    marray = np.ma.masked_where(vlabel - VOIsL != 0, data_mn).mask
    dat1 = np.ma.array(data_mn, mask = marray).compressed()
    marray = np.ma.masked_where(vlabel - VOIsR != 0, data_mn).mask
    dat2 = np.ma.array(data_mn, mask = marray).compressed()

    array_mn = 0.5 * (dat1 + dat2)

    marray = np.ma.masked_where(vlabel - VOIsL != 0, data_sd).mask
    dat3 = np.ma.array(data_sd, mask = marray).compressed()
    marray = np.ma.masked_where(vlabel - VOIsR != 0, data_sd).mask
    dat4 = np.ma.array(data_sd, mask = marray).compressed()

    array_sd = np.sqrt( 0.5 * (dat3 * dat3 + dat4 * dat4) + 0.25 * (dat1 - dat2) * (dat1 - dat2) )
    
    return array_mn, array_sd

def main():
    
    img = nib.load('../../dump.0.nii')
    data = img.get_fdata()
    
    Nx[0], Nx[1], Nx[2], dum, Nx[3] = data.shape
    dx[0], dx[1], dx[2], dt, dum = img.header.get_zooms()
    
    for i in range(3):
        dx[i] *= scale_fac
        lx[i] = Nx[i] * dx[i]

    # dim1
    Npoints = Nx[0] * Nx[2]
    x.append(np.zeros(Npoints, dtype=np.float32))
    y.append(np.zeros(Npoints, dtype=np.float32))

    # dim 1
    c = 0
    for i in range(Nx[0]):
        dum = dx[0] * i
        for j in range(Nx[2]):
            x[0][c] = dum
            y[0][c] = dx[2] * j
            
            c += 1

    dataf = data[:,:,:,0,5].flatten()
    marray = np.ma.masked_where(dataf != 1, dataf).mask
    dat1 = np.ma.array(dataf, mask = marray).compressed()
    marray = np.ma.masked_where(dataf != 2, dataf).mask
    dat2 = np.ma.array(dataf, mask = marray).compressed()
    
    tot_vol = dx[0] * dx[1] * dx[2] * (len(dat1) + len(dat2))

    # time, voi_label, voi_name, vol, M_mn, M_std,  N_mn, N_std,  S_mn, S_std,  F_mn, F_std,  A_mn, A_std,  TYPE_mn, TYPE_std,  VOI_mn, VOI_std,  atrophy_mn, atrophy_std
    #    0          1         2    3     4      5      6      7      8      9     10     11     12     13        14        15       16       17           18           19
    step = np.genfromtxt(fname, dtype=np.int, comments='#',usecols=(0))
    vlabel = np.genfromtxt(fname, dtype=np.int, comments='#',usecols=(1))
    vname = np.genfromtxt(fname, dtype="|U4", comments='#',usecols=(2))
    vol = np.genfromtxt(fname, dtype=np.float, comments='#',usecols=(3))
    
    agent_mn = []
    agent_sd = []
    
    for i in range(len(ag_sgn)):
        c = 2*i + 4
        agent_mn.append(np.genfromtxt(fname, dtype=np.float, comments='#',usecols=(c)))
        c = 2*i + 5
        agent_sd.append(np.genfromtxt(fname, dtype=np.float, comments='#',usecols=(c)))
    
    ### extract the desired VOIs
    time, dum = extract(step,np.zeros_like(step), vlabel, VOIs[0][1], VOIs[0][2])
    time *= scale_time
    
    voi_vol = np.zeros(len(VOIs),dtype=float)
    voi_max = np.zeros((len(VOIs),len(ag_sgn),2),dtype=float)

    voi_ag_mn = []
    voi_ag_sd = []
    for i in range(len(VOIs)):
        vv, dum = extract(vol,np.zeros_like(vol), vlabel, VOIs[i][1], VOIs[i][2])
        voi_vol[i] = np.mean(vv) / tot_vol

        voi_ag_mn.append([])
        voi_ag_sd.append([])
        for j in range(len(ag_sgn)):
            a_mn, a_sd = extract(agent_mn[j],agent_sd[j], vlabel, VOIs[i][1], VOIs[i][2])
            
            voi_ag_mn[i].append(a_mn)
            voi_ag_sd[i].append(a_sd)
            
            voi_max[i,j,0] = np.max(a_mn)
            voi_max[i,j,1] = a_sd[np.argmax(a_mn)]
            
        voi_id[i] = i
        voi_tit[i] = VOIs[i][0]
        voi_col[i] = VOIs[i][5]
        
    make_plot(time, voi_vol, voi_ag_mn, voi_ag_sd, voi_max, data[:,:,:,0,6])
    
#########################
## plot
#########################
def make_plot(time, voi_vol, voi_ag_mn, voi_ag_sd, voi_max, voi_map):
    
    # Plot properties
    n_X = 3
    n_Y = 2
    ax_len_x = 0.7/n_X         # Length of one subplot square box
    ax_len_y = 0.355/n_Y         # Length of one subplot square box
    ax_bx = 0.08               # Beginning/offset of the subplot in the box
    ax_by = 0.54               # Beginning/offset of the subplot in the box
    ax_sepx = 0.10             # Separation length between two subplots
    ax_sepy = 0.065             # Separation length between two subplots
    total_subplots_in_x = n_X  # Total number of subplots

    ax = []
    fig = pl.figure(figsize=(10,12))
    subp = Subplots(fig, ax_len_x, ax_len_y, ax_sepx, ax_sepy, ax_bx, ax_by, total_subplots_in_x)
    c = 0

    mlx   = MultipleLocator(1)
    mly   = MultipleLocator(0.1)
    
    # Panel numbering properties
    props = dict(boxstyle='round', facecolor='white', edgecolor='none', alpha=1.0)

    # empty
    ax.append( subp.addSubplot() )
    ax_id = 0

    ax[ax_id].set_axis_off()

    # fAb
    me = 3
    ax.append( subp.addSubplot() )
    ax_id = 1
    
    val_max = np.max(voi_max[:,me,0])
    for i in range(len(VOIs)):
        ax[ax_id].plot(time, voi_ag_mn[i][me] / val_max,c=VOIs[i][5],lw=VOIs[i][6], ls=VOIs[i][7], label = VOIs[i][0])

    ax[ax_id].set_xlabel(r'$\rm time ~(yr)$')
    ax[ax_id].xaxis.set_minor_locator(mlx)
    ax[ax_id].set_title(ag_tit[me],fontsize=20)
    
    ax[ax_id].set_xlim(0, 13.9)

    # ast
    me = 4
    ax.append( subp.addSubplot() )
    ax_id = 2
    
    for i in range(len(VOIs)):
        ax[ax_id].plot(time, voi_ag_mn[i][me],c=VOIs[i][5],lw=VOIs[i][6], ls=VOIs[i][7], label = VOIs[i][0])

    ax[ax_id].set_xlabel(r'$\rm time ~(yr)$')
    ax[ax_id].xaxis.set_minor_locator(mlx)
    ax[ax_id].set_title(ag_tit[me],fontsize=24)
    
    ax[ax_id].set_xlim(0, 13.9)
    ax[ax_id].set_ylim(0,1)

    # neu
    me = 1
    ax.append( subp.addSubplot() )
    ax_id = 3
    
    val_max = np.max(voi_max[:,me,0])
    for i in range(len(VOIs)):
        ax[ax_id].plot(time, voi_ag_mn[i][me] / val_max,c=VOIs[i][5],lw=VOIs[i][6], ls=VOIs[i][7], label = VOIs[i][0])

    ax[ax_id].xaxis.set_minor_locator(mlx)
    ax[ax_id].set_title(ag_tit[me],fontsize=20)
    
    ax[ax_id].set_xlim(0, 13.9)

    # atrophy
    me = 7
    ax.append( subp.addSubplot() )
    ax_id = 4
    
    for i in range(len(VOIs)):
        ax[ax_id].plot(time, voi_ag_mn[i][me],c=VOIs[i][5],lw=VOIs[i][6], ls=VOIs[i][7], label = VOIs[i][0])

    ax[ax_id].xaxis.set_major_formatter(pl.NullFormatter())
    ax[ax_id].xaxis.set_minor_locator(mlx)
    ax[ax_id].set_title(ag_tit[me],fontsize=20)
    
    ax[ax_id].set_xlim(0, 13.9)
    ax[ax_id].set_ylim(0, 1)
    
    # mic
    me = 0
    ax.append( subp.addSubplot() )
    ax_id = 5
    
    for i in range(len(VOIs)):
        ax[ax_id].plot(time, voi_ag_mn[i][me] / voi_ag_mn[i][me][0],c=VOIs[i][5],lw=VOIs[i][6], ls=VOIs[i][7], label = VOIs[i][0])

    ax[ax_id].xaxis.set_major_formatter(pl.NullFormatter())
    ax[ax_id].xaxis.set_minor_locator(mlx)
    ax[ax_id].set_title(ag_tit[me],fontsize=20)
    
    ax[ax_id].set_xlim(0, 13.9)

    # voi volumes
    ax2 = fig.add_axes([0.1, 0.075, 0.89, 0.12], facecolor=(1.,1.,1.))
    ax2.bar(voi_id, voi_vol+0.01, width=0.8,bottom=-0.01, align='center', color = voi_col, linewidth=0.5, edgecolor='k')
    ax2.set_xticks(np.arange(len(VOIs)))
    ax2.set_xticklabels(voi_tit, rotation=45, ha='center', va='top')
    ax2.axhline(y=0,linewidth=0.5, color='k')
    ax2.set_ylabel(r'$V/V_{\rm par}$',fontsize=20)

    # voi Fmax
    ax3 = fig.add_axes([0.1, 0.205, 0.89, 0.12], facecolor=(1.,1.,1.))
    ax3.bar(voi_id, voi_max[:,3,0], yerr=[np.zeros_like(voi_max[:,3,1]), voi_max[:,3,1]],
            width=0.8,bottom=0, align='center', color = voi_col, linewidth=0.5, edgecolor='k', error_kw=dict(ecolor='k', elinewidth=0.5, capthick=0.5, capsize=3))
    ax3.set_xticks(np.arange(len(VOIs)))
    pl.setp(ax3.get_xticklabels(), visible=False)
    ax3.set_ylabel(r'$F_{\rm max} ~{\rm (\mu M)}$',fontsize=24)

    # voi neu_0
    ax4 = fig.add_axes([0.1, 0.335, 0.89, 0.12], facecolor=(1.,1.,1.))
    ax4.bar(voi_id, voi_max[:,1,0] / np.max(voi_max[:,1,0]), yerr=[np.zeros_like(voi_max[:,1,1]), voi_max[:,1,1] / np.max(voi_max[:,1,0]) ],
            width=0.8,bottom=0, align='center', color = voi_col, linewidth=0.5, edgecolor='k', error_kw=dict(ecolor='k', elinewidth=0.5, capthick=0.5, capsize=3))
    ax4.set_xticks(np.arange(len(VOIs)))
    pl.setp(ax4.get_xticklabels(), visible=False)
    ax4.set_ylim(0,1.7)
    ax4.set_ylabel(r'$N_{\rm o}/N_{\rm max}$',fontsize=20)

    # brain section
    ax5 = fig.add_axes([0.04, 0.48, 0.30, 0.25], facecolor=(1.,1.,1.))
    ax5.set_axis_off()

    import matplotlib.image as mpimg
    im = mpimg.imread('group.png')
    ax5.imshow(im)

    #labeling
    fig.text(0.01,0.975,r'a',fontsize=26)
    fig.text(0.35,0.975,r'b',fontsize=26)
    fig.text(0.67,0.975,r'c',fontsize=26)
    
    fig.text(0.01,0.73,r'd',fontsize=26)
    fig.text(0.35,0.73,r'e',fontsize=26)
    fig.text(0.67,0.73,r'f',fontsize=26)

    fig.text(0.11,0.425,r'g',fontsize=26)
    fig.text(0.11,0.3,r'h',fontsize=26)
    fig.text(0.11,0.17,r'i',fontsize=26)


    """
    ## dim1
    bx = 0.0
    by = 0.5
    xls = np.linspace(bx,bx+lx[0],Nx[0])
    yls = np.linspace(by,by+lx[2],Nx[2])
    x_gr, y_gr = np.meshgrid(xls, yls)
    
    sec = 60
    z = voi_map[:,sec,:].flatten()
    
    # grid the data.
    z_gr = griddata((x[0]+bx,y[0]+by), z, (x_gr, y_gr), method='nearest')

    cmap = matplotlib.colors.ListedColormap(voi_col)

    contour_levels = len(VOIs)
    cs = ax5.contourf(x_gr,y_gr,z_gr, contour_levels, cmap=cmap) # levels = np.linspace(cut_min,cut_max,contour_levels), extend='both') # , norm = LogNorm())
    
    cs.cmap.set_under('None')
    cs.cmap.set_over('None')
    """

    pl.savefig("voi_time.png", format = 'png', dpi=200, orientation='landscape')
    pl.show()
    #sys.stdin.read(1)
    pl.close()


if __name__ == '__main__':
    main()
