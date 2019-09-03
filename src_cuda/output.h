#ifndef OUTPUT_H
#define OUTPUT_H

#include "pointers.h"

using namespace std;

class Output {
 public:
  Output();
  ~Output();

  void lammpstrj(class VirtualBrain*);
  void restart(class VirtualBrain*);
  void statistics(class VirtualBrain*);
  void statistics_sphere(class VirtualBrain*);
  void dump(class VirtualBrain*);

  void dump_txt(VirtualBrain*, vector<string>);
  void dump_mri(VirtualBrain*, vector<string>);
  void dump_mricrogl(VirtualBrain*, vector<string>);

  void sort_tag(VirtualBrain*, double*, int);
  tagint find_tag(VirtualBrain*, double, double, double);
  //nifti_image *nifti_image_setup(VirtualBrain*, vector<string>, int);
  nifti_image *nifti_image_setup(VirtualBrain*, vector<string>, const int[], int);

  fstream fb;
  FILE* fw;
  string rname, sname;
  bool do_dump, do_statistics;
  int revery,severy;

  vector<vector<string>> dump_arg;

};

#endif
