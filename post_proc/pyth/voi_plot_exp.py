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
VOIs.append( (r'$\rm FrG$', 210, 211,  'FLG ', 'FRG ', 'blue',          2.2, '--' ) )
VOIs.append( (r'$\rm TeG$', 218, 219,  'TLG ', 'TRG ', 'green',         2.0, '--' ) )
VOIs.append( (r'$\rm PaG$', 6, 2,      'PLG ', 'PRG ', 'magenta',       1.9, '--' ) )
VOIs.append( (r'$\rm OcG$', 8, 4,      'OLG ', 'ORG ', 'lightgreen',    1.8, '--' ) )
VOIs.append( (r'$\rm Cau$', 39, 53,    'CaL ', 'CaR ', 'orange',        1.7, '--' ) )
VOIs.append( (r'$\rm Put$', 14, 16,    'PuL ', 'PuR ', 'steelblue',     1.6, '-' ) )
VOIs.append( (r'$\rm Tha$', 102, 203,  'ThL ', 'ThR ', 'hotpink',       1.6, '--' ) )
VOIs.append( (r'$\rm Cer$', 67, 76,    'CeL ', 'CeR ', 'firebrick',     1.7, '-' ) )

# Fig 2: arnold1991topographical
# Neurofibrillary tangle (NFT) values
NFTexpO = (r'$\rm OcX_a$', 1.25, 2.37)
NFTexpP = (r'$\rm PaX_a$', 1.66, 2.90)
NFTexpF = (r'$\rm FrX_a$', 2.01, 3.19)
NFTexpT = (r'$\rm TeX_a$', 2.64, 3.77)
NFTexpL = (r'$\rm LiX_a$', 2.92, 3.62)
# Neuritic plaque (NP) values
NPexpO = (r'$\rm OcX_a$', 2.62, 3.38)
NPexpP = (r'$\rm PaX_a$', 2.38, 3.16)
NPexpF = (r'$\rm FrX_a$', 2.19, 2.66)
NPexpT = (r'$\rm TeX_a$', 2.56, 3.14)
NPexpL = (r'$\rm LiX_a$', 1.85, 2.26)

# Fig 3: mclean1999soluble
# Insoluble Amyloid beta (ug/g)
NPexp2O = (r'$\rm OcX_m$', 22.84, 48.32)
NPexp2F = (r'$\rm FrX_m$', 21.18, 32.10)
NPexp2T = (r'$\rm TeX_m$', 29.62, 68.18)
NPexp2I = (r'$\rm InX_m$', 9.10, 17.21)
NPexp2Hip = (r'$\rm HipX_m$', 12.41, 32.43)
NPexp2Put = (r'$\rm PutX_m$', 7.61, 16.38)
NPexp2Tha = (r'$\rm ThaX_m$', 8.11, 20.85)
NPexp2Cer = (r'$\rm CerX_m$', 11.42, 30.45)

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
    
    x = []
    y = []

    Nx = np.zeros(4, dtype=np.int32)
    dx = np.zeros(3, dtype=np.float32)
    lx = np.zeros(3, dtype=np.float32)
    
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
            
            voi_max[i,j,0] = a_mn[50] #np.max(a_mn)
            voi_max[i,j,1] = a_sd[50] #a_sd[np.argmax(a_mn)]
            
        voi_id[i] = i
        voi_tit[i] = VOIs[i][0]
        voi_col[i] = VOIs[i][5]

        print(voi_id[i], voi_tit[i], voi_vol[i], voi_max[i,3,0])

    # figure preparation
    bw = 1.0
    sepx = 1.0
    dx = 0.02

    lbl = () #np.zeros(len(VOIs) + len(NPexp), dtype='|U10')
    y   = np.zeros(19, dtype=float)
    yerr= np.zeros(19, dtype=float)
    x   = np.zeros(19, dtype=float)

    rgb = np.zeros((19, 3), dtype=float)
    
    norm = voi_max[3,3,0]
    normex = NPexpO[1]
    normex2 = NPexp2O[1]
    
    i = 0
    # Temporal lobe
    y[i+0] = voi_max[1,3,0] / norm
    y[i+1] = NPexpT[1] / normex
    y[i+2] = NPexp2T[1] / normex2
    
    yerr[i+0] = voi_max[1,3,1] / norm
    yerr[i+1] = (NPexpT[2] - NPexpT[1]) / normex
    yerr[i+2] = (NPexp2T[2] - NPexp2T[1]) / normex2

    lbl += ( voi_tit[1],)
    lbl += ( NPexpT[0],)
    lbl += ( NPexp2T[0],)

    rgb[i+0,0], rgb[i+0,1], rgb[i+0,2] = 0.5, 0.5, 0.5
    rgb[i+1,0], rgb[i+1,1], rgb[i+1,2] = 0.7, 0.7, 0.7
    rgb[i+2,0], rgb[i+2,1], rgb[i+2,2] = 0.9, 0.9, 0.9
    
    x[i+0] = 0.0
    x[i+1] = x[i+0] + bw + dx
    x[i+2] = x[i+1] + bw + dx

    i += 3

    # Frontal lobe
    y[i+0] = voi_max[0,3,0] / norm
    y[i+1] = NPexpF[1] / normex
    y[i+2] = NPexp2F[1] / normex2
    
    yerr[i+0] = voi_max[0,3,1] / norm
    yerr[i+1] = (NPexpF[2] - NPexpF[1]) / normex
    yerr[i+2] = (NPexp2F[2] - NPexp2F[1]) / normex2

    lbl += ( voi_tit[0],)
    lbl += ( NPexpF[0],)
    lbl += ( NPexp2F[0],)
    
    x[i+0] = x[i-1] + sepx + bw
    x[i+1] = x[i+0] + bw + dx
    x[i+2] = x[i+1] + bw + dx

    rgb[i+0,0], rgb[i+0,1], rgb[i+0,2] = 0.5, 0, 0
    rgb[i+1,0], rgb[i+1,1], rgb[i+1,2] = 0.7, 0, 0
    rgb[i+2,0], rgb[i+2,1], rgb[i+2,2] = 0.9, 0, 0

    i += 3

    # Occipital lobe
    y[i+0] = voi_max[3,3,0] / norm
    y[i+1] = NPexpO[1] / normex
    y[i+2] = NPexp2O[1] / normex2
    
    yerr[i+0] = voi_max[3,3,1] / norm
    yerr[i+1] = (NPexpO[2] - NPexpO[1]) / normex
    yerr[i+2] = (NPexp2O[2] - NPexp2O[1]) / normex2

    lbl += ( voi_tit[3],)
    lbl += ( NPexpO[0],)
    lbl += ( NPexp2O[0],)
    
    x[i+0] = x[i-1] + sepx + bw
    x[i+1] = x[i+0] + bw + dx
    x[i+2] = x[i+1] + bw + dx
    
    rgb[i+0,0], rgb[i+0,1], rgb[i+0,2] = 0, 0.5, 0
    rgb[i+1,0], rgb[i+1,1], rgb[i+1,2] = 0, 0.7, 0
    rgb[i+2,0], rgb[i+2,1], rgb[i+2,2] = 0, 0.9, 0
    
    i += 3

    # Parietal lobe
    y[i+0] = voi_max[2,3,0] / norm
    y[i+1] = NPexpP[1] / normex
    
    yerr[i+0] = voi_max[2,3,1] / norm
    yerr[i+1] = (NPexpP[2] - NPexpP[1]) / normex

    lbl += ( voi_tit[2],)
    lbl += ( NPexpP[0],)

    x[i+0] = x[i-1] + sepx + bw
    x[i+1] = x[i+0] + bw + dx
    
    rgb[i+0,0], rgb[i+0,1], rgb[i+0,2] = 0, 0, 0.5
    rgb[i+1,0], rgb[i+1,1], rgb[i+1,2] = 0, 0, 0.7
    
    i += 2

    # inners
    y[i+0] = voi_max[5,3,0] / norm
    y[i+1] = voi_max[6,3,0] / norm
    y[i+2] = voi_max[7,3,0] / norm
    
    y[i+3] = NPexp2Put[1] / normex2
    y[i+4] = NPexp2Tha[1] / normex2
    y[i+5] = NPexp2Cer[1] / normex2
    
    yerr[i+0] = voi_max[5,3,1] / norm
    yerr[i+1] = voi_max[6,3,1] / norm
    yerr[i+2] = voi_max[7,3,1] / norm
    
    yerr[i+3] = (NPexp2Put[2] - NPexp2Put[1]) / normex2
    yerr[i+4] = (NPexp2Tha[2] - NPexp2Tha[1]) / normex2
    yerr[i+5] = (NPexp2Cer[2] - NPexp2Cer[1]) / normex2
    
    lbl += ( voi_tit[5], )
    lbl += ( voi_tit[6], )
    lbl += ( voi_tit[7], )
    
    lbl += ( NPexp2Put[0], )
    lbl += ( NPexp2Tha[0], )
    lbl += ( NPexp2Cer[0], )

    x[i+0] = x[i-1] + sepx + bw
    x[i+1] = x[i+0] + bw + dx
    x[i+2] = x[i+1] + bw + dx
    
    x[i+3] = x[i+2] + bw + dx + 5.*dx
    x[i+4] = x[i+3] + bw + dx
    x[i+5] = x[i+4] + bw + dx
    
    rgb[i+0,0], rgb[i+0,1], rgb[i+0,2] = 0.5, 0, 0.5
    rgb[i+1,0], rgb[i+1,1], rgb[i+1,2] = 0.6, 0, 0.6
    rgb[i+2,0], rgb[i+2,1], rgb[i+2,2] = 0.7, 0, 0.7
    rgb[i+3,0], rgb[i+3,1], rgb[i+3,2] = 0.8, 0, 0.8
    rgb[i+4,0], rgb[i+4,1], rgb[i+4,2] = 0.9, 0, 0.9
    rgb[i+5,0], rgb[i+5,1], rgb[i+5,2] = 1, 0, 1
    
    i += 6

    # Limbic lobe / Insula cortex
    y[i+0] = NPexpL[1] / normex
    y[i+1] = NPexp2I[1] / normex2
    
    yerr[i+0] = (NPexpL[2] - NPexpL[1]) / normex
    yerr[i+1] = (NPexp2I[2] - NPexp2I[1]) / normex2
    
    lbl += ( NPexpL[0], )
    lbl += ( NPexp2I[0], )

    x[i+0] = x[i-1] + sepx + bw
    x[i+1] = x[i+0] + bw + dx

    rgb[i+0,0], rgb[i+0,1], rgb[i+0,2] = 0, 0.5, 0.5
    rgb[i+1,0], rgb[i+1,1], rgb[i+1,2] = 0, 0.7, 0.7
    
    i += 2
    
    make_plot(x, y, yerr, lbl, rgb, bw)
    
#########################
## plot
#########################
def make_plot(x, y, yerr, lbl, rgb,bw):
    fig = pl.figure(figsize=(10,4))
    
    ax = fig.add_axes([0.1, 0.18, 0.89, 0.8], facecolor=(1.,1.,1.))
    
    ax.bar(x, y, yerr= [np.zeros_like(yerr), yerr],
           width=bw, bottom=0, align='center', color = rgb, linewidth=0.5, edgecolor='k', error_kw=dict(ecolor='k', elinewidth=0.5, capthick=0.5, capsize=3))

    print(lbl)
    ax.set_xticks(x)
    ax.set_xticklabels(lbl, rotation=45, ha='center', va='top')
    ax.set_ylim(0,2)
    ax.set_ylabel(r'$\rm A\beta ~plaques$',fontsize=20)

    pl.savefig("NP_exp.png", format = 'png', dpi=200, orientation='landscape')
    pl.show()
    #sys.stdin.read(1)
    pl.close()

if __name__ == '__main__':
    main()
