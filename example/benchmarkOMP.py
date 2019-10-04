###########
import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import numpy as np
import matplotlib.pyplot as pl
from scipy.optimize import fsolve, curve_fit
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import axes3d
from math import sqrt, floor, pi
from matplotlib import rc, rcParams, cm
from matplotlib.ticker import MultipleLocator
import matplotlib.mlab as mlab
from matplotlib.mlab import griddata
import subprocess

from tempfile import mkstemp
from shutil import move
from os import fdopen, remove

import sys

rc('font',**{'family':'sans-serif','serif':['Times'], 'weight': 'bold', 'size' : 20})
rc('text', usetex=True)
rcParams.update({'figure.autolayout': True})

def replace(file_path, line_tit, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                ls = line.split()
                if (len(ls) > 0) :
                    if (ls[0] == 'partition'):
                        line = subst
                new_file.write(line)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)
#Nproc = ( ( 1, 1, 1, 1))
Nproc = ( ( 1, 1,  1, 1),
          ( 2, 1,  2, 1),
          ( 4, 1,  4, 1),
          ( 6, 1,  6, 1),
          ( 8, 1,  8, 1),
          (10, 1, 10, 1),
          (12, 2,  6, 1),
          ( 2, 2,  1, 1),
          ( 4, 2,  2, 1),
          ( 8, 2,  4, 1) )

#Nthread = (1, 2, 4, 8, 12, 16, 20, 24)
Nthread =2
method = ('OMP', )
exes = ('../src_cuda/braimmu.exe', )


def main():
    

    mu_speed = np.zeros((len(Nproc),len(method)),dtype=np.float32)
    sig_speed = np.zeros((len(Nproc),len(method)),dtype=np.float32)

    fwname = 'benchmarkOMP+MPI.out'
    fw = open(fwname, 'w')
    fw.write('# benchmark results \n')
    line = 'omp_threads ncore nx ny nz method mu sig\n'
    fw.write(line)
    fw.close()
    
   
    for sim in range(len(method)):
        for it in range(len(Nproc)):
           # for j in range (len(Nthread)):
               ncore, nx, ny, nz = Nproc[it]
            
               subst = 'partition %i %i %i \n' % (nx,ny,nz)
            
            # setup the partitions in the input file
               fname = 'connOMP.run'
               replace(fname,'partition', subst)

               print(ncore, nx,ny,nz, method[sim])
            
            #command = 'mpirun -use-hwthread-cpus -np %i %s %s' % (ncore, exes[sim], fname)
               command = 'export OMP_NUM_THREADS=%i; srun -u -n %i %s connectome %s' % (Nthread, ncore, exes[sim], fname)
               print(command)
               p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
               p.wait()
            
               fr = open ( 'log.braimmu',"r")
               llist = fr.readlines()
               fr.close()

               print(llist[-1])
            
               ls = llist[-1].split()
               mu_speed[it,sim] = ls[2]
               sig_speed[it,sim] = ls[4]
            
               print(mu_speed[it,sim], sig_speed[it,sim])
            
               fw = open(fwname, 'a')
               line = '%i %i %i %i %i %s %g %g\n' % (Nthread,ncore,nx,ny,nz,method[sim],mu_speed[it,sim],sig_speed[it,sim])
               fw.write(line)
               fw.close()

if __name__ == '__main__':
    main()
