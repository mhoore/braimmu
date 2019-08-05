#include "math.h"
#include "cstring"

#include "iomanip"
#include "iostream"
#include "fstream"
#include "stdio.h"
#include "stdlib.h"
#include "string"
#include "sstream"

#include "nifti1.h"
#include "nifti1_io.h"
#include "znzlib.h"

int main(int narg, char **arg) {

  if (narg < 2) {
    printf ("ERROR: arg1 = nifti image file \n");
    exit (EXIT_FAILURE);
  }

  nifti_image *nim = nifti_image_read(arg[1],1);

  if(!nim) {
    printf("Error: cannot open nifti image %s \n", arg[1]);
    exit(1);
  }

  printf("##################### \n");
  printf("NIFTI image %s is read. \n", arg[1]);
  printf("NIFTI image properties: ");
  printf("ndim = %i \n", nim->ndim);
  for (int i=1; i<8; i++)
    printf("dim[%i] = %i, pixdim[i] = %g \n",
           i,nim->dim[i],i,nim->pixdim[i]);
  printf("nvox = %lli \n", nim->nvox);
  printf("nbyper = %i \n", nim->nbyper);
  printf("datatype = %i \n", nim->datatype);

  printf("calmin = %g, calmax = %g \n", nim->cal_min, nim->cal_max);
  printf("toffset = %g \n", nim->toffset);

  printf("xyz_units = %i, time_units = %i \n", nim->xyz_units, nim->time_units);
  printf("nifti_type = %i \n", nim->nifti_type);

  printf("intent_code = %i \n", nim->intent_code);

  printf("description: %s \n", nim->descrip);
  printf("##################### \n");

  return 0;
}
