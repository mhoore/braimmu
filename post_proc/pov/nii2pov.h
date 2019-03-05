#ifndef nii_to_pov_H
#define nii_to_pov_H

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

using namespace std;

/// in case number of voxels is more than 2 bilions
typedef int64_t tagint;
#ifndef PRId64
#define PRId64 "ld"
#endif
#define TAGINT_FORMAT "%" PRId64
/// otherwise
//#define TAGINT_FORMAT "%d"
//typedef int32_t tagint;

class N2P {
 public:

  N2P(int, char **);
  ~N2P();

  void *smalloc(size_t n, const char *);
  void sfree(void *);

  void read_image(char*);
  int mri_boundaries(nifti_image*);
  void mri_topology(nifti_image*);
  void settings ();
  void file_write ();

  double interpolate(double, double, double, double, double);
  void jet(double, double*);

  nifti_image *nim;

  double vlen[3], lbox[3], boxlo[3], boxhi[3];
  int nv[3];

  double minval, maxval; // minimum and maximum values to be shown
  int sec_id; // which section to show

  uint8_t * ptr8;
  int16_t *ptr16;
  int32_t *ptr32;
  float *ptrf;
  double *ptrd;

  double **data, **x;
  int dsize, nvoxel;

  vector<string> darg;

  int me;

  double cam[3], look[3], right[3], up[3];

  char fname[FILENAME_MAX];
  FILE* fw;


  /* ----------------------------------------------------------------------
      create a 1d array
  ------------------------------------------------------------------------- */
    template <typename TYPE>
    TYPE *create(TYPE *&array, int n, const char *name) {
      tagint nbytes = ((tagint) sizeof(TYPE)) * n;
      array = (TYPE *) smalloc(nbytes,name);
      return array;
    }

  /* ----------------------------------------------------------------------
     destroy a 1d array
  ------------------------------------------------------------------------- */
    template <typename TYPE>
    void destroy(TYPE *&array) {
      sfree(array);
      array = NULL;
    }

  /* ----------------------------------------------------------------------
     create a 2d array
  ------------------------------------------------------------------------- */
    template <typename TYPE>
    TYPE **create(TYPE **&array, int n1, int n2, const char *name) {
      tagint nbytes = ((tagint) sizeof(TYPE)) * n1*n2;
      TYPE *data = (TYPE *) smalloc(nbytes,name);
      nbytes = ((tagint) sizeof(TYPE *)) * n1;
      array = (TYPE **) smalloc(nbytes,name);

      tagint n = 0;
      for (int i = 0; i < n1; i++) {
        array[i] = &data[n];
        n += n2;
      }
      return array;
    }

  /* ----------------------------------------------------------------------
     destroy a 2d array
  ------------------------------------------------------------------------- */
    template <typename TYPE>
    void destroy(TYPE **&array) {
      if (array == NULL) return;
      sfree(array[0]);
      sfree(array);
      array = NULL;
    }

  /* ----------------------------------------------------------------------
     destroy a 2d array with 2nd index offset
  ------------------------------------------------------------------------- */
    template <typename TYPE>
    void destroy2d_offset(TYPE **&array, int offset) {
      if (array == NULL) return;
      sfree(&array[0][offset]);
      sfree(array);
      array = NULL;
    }

  /* ----------------------------------------------------------------------
     create a 3d array
  ------------------------------------------------------------------------- */
    template <typename TYPE>
    TYPE ***create(TYPE ***&array, int n1, int n2, int n3, const char *name) {
      tagint nbytes = ((tagint) sizeof(TYPE)) * n1*n2*n3;
      TYPE *data = (TYPE *) smalloc(nbytes,name);
      nbytes = ((tagint) sizeof(TYPE *)) * n1*n2;
      TYPE **plane = (TYPE **) smalloc(nbytes,name);
      nbytes = ((tagint) sizeof(TYPE **)) * n1;
      array = (TYPE ***) smalloc(nbytes,name);

      int i,j;
      tagint m;
      tagint n = 0;
      for (i = 0; i < n1; i++) {
        m = ((tagint) i) * n2;
        array[i] = &plane[m];
        for (j = 0; j < n2; j++) {
          plane[m+j] = &data[n];
          n += n3;
        }
      }
      return array;
    }

  /* ----------------------------------------------------------------------
     destroy a 3d array
  ------------------------------------------------------------------------- */
    template <typename TYPE>
    void destroy(TYPE ***&array) {
      if (array == NULL) return;
      sfree(array[0][0]);
      sfree(array[0]);
      sfree(array);
      array = NULL;
    }

};

#endif
