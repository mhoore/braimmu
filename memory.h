#ifndef MEMORY_H
#define MEMORY_H

#include <vector>
#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#include "pointers.h"

namespace brain_NS {

class Memory {
 public:
  Memory();
  ~Memory();

  // union data struct for packing 32-bit and 64-bit ints into double bufs
  // this avoids aliasing issues by having 2 pointers (double,int)
  union ubuf {
    double d;
    int64_t i;
    ubuf(double arg) : d(arg) {}
    ubuf(int64_t arg) : i(arg) {}
    ubuf(int arg) : i(arg) {}
  };

  void *smalloc(size_t n, const char *);
  void sfree(void *);

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

}

#endif
