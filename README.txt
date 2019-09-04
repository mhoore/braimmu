Here you see how you can run a simulation:

For running this code, your system needs message passing interface (mpi) compiler on c++ (mpicxx) to be already installed on your machine.

1- Go to ./src folder and enter the below command
>> make

2- If the code compiled well, a file named "braimmu.exe" would be generated.

3- Run an example simulation using one of input scripts, in.run or in.cuda in ./example/ directory:
Note: number of cores should match the partitions specifications inside the input script
>> mpirun -np 4 ../src/braimmu.exe connectome conn.run

4- The brain atlases are used as initial conditions, and taken from multiple sources:

* Neuroimaging and Surgical Technologies Lab, "http://nist.mni.mcgill.ca/?p=904":
  ICBM 2009c Nonlinear Symmetric 1×1x1mm template: NIFTI 55MB
  ./base/wm.nii  --> white matter: mni_icbm152_wm_tal_nlin_sym_09c.nii
  ./base/gm.nii  --> gray matter: mni_icbm152_gm_tal_nlin_sym_09c.nii
  ./base/csf.nii --> cerebrospinal fluid: mni_icbm152_csf_tal_nlin_sym_09c.nii

* Courtesy of Dr. Vladimir S. Fonov, McGill University, Montreal Neurological Institute, Brain Imaging Center:
  ./base/voi.nii --> volume of interest (VOI): mni_icbm152_lob_tal_nlin_sym_09c.nii
  ./base/voi.txt --> VOIs descriptions

* UCLA Brain Mapping Center, "http://www.bmap.ucla.edu/portfolio/atlases/ICBM_DTI-81_Atlas/":
  ICBM DTI-81 Atlas
  ./base/rgb.nii --> neural connectome: ICBM_COLOR.nii

5- Enjoy!
