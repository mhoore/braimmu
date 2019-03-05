Here you see how you can run a simulation:

For running this code, your system needs message passing interface (mpi) compiler on c++ (mpicxx) to be already installed on your machine.

1- Go to src folder and enter the below command
>> make

2- If the code compiled well, a file named "braimmu.exe" would be generated.

3- Enter the below command to run the code on 12 cores (You can set the number of cores in the file "in.run"):
>> mpirun -use-hwthread-cpus -np 12 ../src/braimmu.exe in.run

4- The initial file is taken from "https://github.com/neurolabusc/MRIcroGL/blob/master/Distro/mni152_2009bet.nii.gz".

5- Enjoy!
