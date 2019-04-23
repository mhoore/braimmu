Here you see how you can run a simulation:

For running this code, your system needs message passing interface (mpi) compiler on c++ (mpicxx) to be already installed on your machine.

1- Go to ./src folder and enter the below command
>> make

2- If the code compiled well, a file named "braimmu.exe" would be generated.

3- Run an example simulation using one of input scripts, in.run or in.cuda in ./example/ directory:
Note: number of cores should match the partitions specifications inside the input script
>> mpirun -use-hwthread-cpus -np 4 ../src/braimmu.exe in.run

4- The nifti images as initial conditions have been taken from "Neuroimaging and Surgical Technologies Lab", with web address:"http://nist.mni.mcgill.ca/?p=904"; the original file names are:
ICBM 2009c Nonlinear Symmetric 1×1x1mm template: NIFTI 55MB
./example/wm.nii  -> mni_icbm152_wm_tal_nlin_sym_09c.nii
./example/gm.nii  -> mni_icbm152_gm_tal_nlin_sym_09c.nii
./example/csf.nii -> mni_icbm152_csf_tal_nlin_sym_09c.nii

5- Enjoy!
