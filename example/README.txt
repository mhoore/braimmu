Here you see how you can run a simulation:

For running this code, your system needs message passing interface (mpi) compiler on c++ (mpicxx) to be already installed on your machine.

1- Go to src folder and enter the below command
>> make

2- If the code compiled well, a file named "braimmu.exe" would be generated.

3- Enter the below command to run the code on 12 cores (You can set the number of cores in the file "in.run"):
>> mpirun -use-hwthread-cpus -np 12 ../src/braimmu.exe in.run

4- The nifti images as initial conditions have been taken from "Neuroimaging and Surgical Technologies Lab", with web address:"http://nist.mni.mcgill.ca/?p=904"; the original file names are:
ICBM 2009c Nonlinear Symmetric 1×1x1mm template: NIFTI 55MB
wm.nii  -> mni_icbm152_wm_tal_nlin_sym_09c.nii
gm.nii  -> mni_icbm152_gm_tal_nlin_sym_09c.nii
csf.nii -> mni_icbm152_csf_tal_nlin_sym_09c.nii

5- Enjoy!
