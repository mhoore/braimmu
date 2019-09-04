#ifndef INIT_H
#define INIT_H

#include "pointers.h"

using namespace std;

class Init {
 public:
  Init();
  ~Init();

  void setup(class VirtualBrain*);
  void read_mri(class VirtualBrain*);
  void print_mri_properties(class VirtualBrain*, nifti_image*, string);
  void boundaries(class VirtualBrain*);
  void voxels(class VirtualBrain*, bool);
  void neighbor(class VirtualBrain*);

  int mri_boundaries(class VirtualBrain*, nifti_image*);

  vector<vector<string>> mri_arg;

  tagint find_tag(class VirtualBrain*, int, int, int);

  int map(tagint itag) {
    return maparr[itag];
  }

private:
  vector<int> maparr; // map from voxel global id to local id

};

#endif
