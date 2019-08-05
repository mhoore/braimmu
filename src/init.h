#ifndef INIT_H
#define INIT_H

#include <vector>
#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#include "nifti1.h"
#include "nifti1_io.h"
#include "znzlib.h"

#include "pointers.h"
#include "brain.h"

namespace brain_NS {

class Init : public Memory {
 public:
  Init();
  ~Init();

  void setup(class Brain*);
  void read_mri(class Brain*);
  void print_mri_properties(class Brain*, nifti_image*, string);
  void boundaries(class Brain*);
  void voxels(class Brain*, int);
  void neighbor(class Brain*);
  void allocations(class Brain*, int);
  void set_parameters(class Brain*);

  int mri_boundaries(class Brain*, nifti_image*);
  void mri_topology(class Brain*, nifti_image*);

  //int map(class Brain*, tagint);

  vector<vector<string>> mri_arg;

private:
  tagint find_tag(class Brain*, int, int, int);

  uint8_t * ptr8;
  int16_t *ptr16;
  int32_t *ptr32;
  float *ptrf;
  double *ptrd;
  uint8_t *ptr_rgb;

  int *map; // map from voxel global id to local id

};

}

#endif
