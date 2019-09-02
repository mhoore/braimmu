#ifndef INIT_H
#define INIT_H

#include "pointers.h"

using namespace std;
using namespace ns_connectome;

class Init {
 public:
  Init();
  ~Init();

  void setup(class VirtualBrain*);
  void read_mri(class VirtualBrain*);
  void print_mri_properties(class VirtualBrain*, nifti_image*, string);
  void boundaries(class VirtualBrain*);
  void voxels(class VirtualBrain*, int);
  void neighbor(class VirtualBrain*);
  void allocations(class VirtualBrain*, int);
  void set_parameters(class VirtualBrain*);

  int mri_boundaries(class VirtualBrain*, nifti_image*);
  void mri_topology(class VirtualBrain*, nifti_image*);

  vector<vector<string>> mri_arg;

private:
  tagint find_tag(class VirtualBrain*, int, int, int);

  uint8_t * ptr8;
  int16_t *ptr16;
  int32_t *ptr32;
  float *ptrf;
  double *ptrd;
  uint8_t *ptr_rgb;

  vector<int> map; // map from voxel global id to local id

};

#endif
