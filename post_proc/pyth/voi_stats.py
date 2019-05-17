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

# Data Input
scale_fac = 1.e-3

if len(sys.argv) < 6:
    print(sys.argv)
    sys.exit("Error: chack the arguments")

sec = []
frname = sys.argv[1]
step_init = int(sys.argv[2])
step_freq = int(sys.argv[3])
step_final = int(sys.argv[4])
fwname = sys.argv[5]

AGENTS = ("mic","neu","sAb","fAb","ast","typ","group")
ag_tit = (r'$\rm Microglia$', r'$\rm Neurons$', r'$\rm [sA\beta]$',r'$\rm [fA\beta]$',r'$\rm Astrogliosis$',r'$\rm Tissue$',r'$\rm VOI$',r'$\rm Atrophy$')
ag_names = (r'$M ~{\rm (mL^{-1})}$', r'$N ~{\rm (mL^{-1})}$', r'$S ~{\rm (\mu M)}$',r'$F ~{\rm (\mu M)}$',r'$\rm Astrogliosis$',r'$\rm Tissue$',r'$\rm VOI$',r'$\rm Atrophy$')
ag_sgn = ('M','N','S','F','A','TYPE','VOI','atrophy')

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

Nx = np.zeros(4, dtype=np.int32)
dx = np.zeros(3, dtype=np.float32)
lx = np.zeros(3, dtype=np.float32)

def main():

    fw = open(fwname, 'w')
    fw.write('# voi data ')
    line = '# time, voi_label, voi_name, vol, '
    for i in range(len(ag_sgn)):
        line += '%s_mn, %s_std,  ' % (ag_sgn[i],ag_sgn[i])
    line += '\n'
    fw.write(line)
    fw.close()

    fr = frname + str(step_init) + '.nii'
    img = nib.load(fr)
    data0 = img.get_fdata()
    
    step = step_init
    while (step <= step_final):
        fw = open(fwname, 'a')
        
        fr = frname + str(step) + '.nii'
        
        img = nib.load(fr)
        data = img.get_fdata()

        Nx[0], Nx[1], Nx[2], dum, Nx[3] = data.shape
        dx[0], dx[1], dx[2], dt, dum = img.header.get_zooms()

        for i in range(3):
            dx[i] *= scale_fac
            lx[i] = Nx[i] * dx[i]

        data_tis = data[:,:,:,0,-1]
        for i in range(len(voi)):
            print("step ", step, " voi ", i, "/",len(voi))
            line = '%i %i %s  ' % (step, voi[i].label, voi[i].name)
            
            marray = np.ma.masked_where(data_tis - voi[i].label != 0, data_tis).mask
            for me in range(Nx[3]):
                data_voi = np.ma.array(data[:,:,:,0,me], mask = marray).compressed()
                
                if (me == 0):
                    vol = dx[0] * dx[1] * dx[2] * len(data_voi)
                    line += '%g  ' % vol
                
                mu = np.mean(data_voi)
                sig = np.std(data_voi)
                
                line += '%g %g  ' % (mu, sig)
            
            # atrophy (1-N)
            me = 1
            data_voi0 = np.ma.array(data0[:,:,:,0,me], mask = marray).compressed()
            data_voi  = np.ma.array( data[:,:,:,0,me], mask = marray).compressed()
            data_voi = - data_voi / data_voi0 + 1.0
            data_voi = ma.masked_invalid(data_voi).compressed()
            
            mu = np.mean(data_voi)
            sig = np.std(data_voi)
            
            line += '%g %g  \n' % (mu, sig)
            fw.write(line)
        
        fw.close()
        step += step_freq

if __name__ == '__main__':
    main()
