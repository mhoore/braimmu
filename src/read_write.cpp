/* ----------------------------------------------------------------------
 * A basic example to read/write a nifti dataset (e.g. cp command).
 *
 * compile example (consider -pedantic or -Wall):
 *
 * gcc -o clib_01_read_write clib_01_read_write.c -I../include -L../lib -lniftiio -lznz -lz -lm
 *
 * R Reynolds   14 Apr 2009
 *----------------------------------------------------------------------
 */

#include <vector>
#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#include <stdint.h>

#include <nifti1.h>
#include <nifti1_io.h>

void lammpstrj(nifti_image*);

int main(int argc, char * argv[])
{
   nifti_image *nim = NULL;
   char *fin = NULL, *fout = NULL;

   /* read input dataset, including data */
   nim = nifti_image_read(fin, 1);
   if( !nim ) {
      fprintf(stderr,"** failed to read NIfTI image from '%s'\n", fin);
      return 2;
   }

   /* assign nifti_image fname/iname pair, based on output filename
      (request to 'check' image and 'set_byte_order' here) */
   if( nifti_set_filenames(nim, fout, 1, 1) ) return 1;

   printf("HERE, dim[0] = %i, %i %i %i %lli \n",nim->ndim,nim->dim[1],nim->dim[2],nim->dim[3],nim->nvox);
   printf("HERE, dim[0] = %i, %i %g %g %g \n",nim->dim[0],nim->nbyper,nim->dx,nim->dy,nim->dz);
   printf("HERE, %s bytes per data = %i \n",nifti_datatype_string(nim->dt),nim->nbyper);

   lammpstrj(nim);

   /* if we get here, write the output dataset */
   nifti_image_write( nim );
   
   /* and clean up memory */
   nifti_image_free( nim );

   return 0;
}



/* ----------------------------------------------------------------------*/
void lammpstrj(nifti_image *nim) {

  int nx = nim->dim[1];
  int ny = nim->dim[2];
  int nz = nim->dim[3];
  int nvox = nim->nvox;
  double dx = nim->dx;
  double dy = nim->dy;
  double dz = nim->dz;

  uint8_t *ptr = (uint8_t *) nim->data;

  FILE* fw;
  fw = fopen("init.lammpstrj","w");

  fprintf(fw,"ITEM: TIMESTEP \n");
  fprintf(fw,"%i \n", 0);

  fprintf(fw,"ITEM: NUMBER OF ATOMS \n");
  fprintf(fw,"%i \n", nvox);
  fprintf(fw,"ITEM: BOX BOUNDS pp pp pp \n");
  fprintf(fw,"%g %g \n", 0.0,nx*dx);
  fprintf(fw,"%g %g \n", 0.0,ny*dy);
  fprintf(fw,"%g %g \n", 0.0,nz*dz);

  fprintf(fw,"ITEM: ATOMS id type x y z \n");
  int c = 0;
  for (int i=0; i<nz; i++) {
    double x = dx * i;
    for (int j=0; j<ny; j++) {
      double y = dy * j;
      for (int k=0; k<nx; k++) {
        double z = dz * k;
        fprintf(fw,"%i ", c+1); // tag
        fprintf(fw,"%u ", ptr[c++]); // type
        fprintf(fw,"%g ", x); // x
        fprintf(fw,"%g ", y); // y
        fprintf(fw,"%g \n", z); // z
        printf("%i \n", c);
      }
    }
  }

  fclose(fw);

}
