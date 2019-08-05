#ifndef OUTPUT_H
#define OUTPUT_H

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

using namespace std;
namespace brain_NS {

class Output {
 public:
  Output();
  ~Output();

  void lammpstrj(class Brain*);
  void restart(class Brain*);
  void statistics(class Brain*);
  void statistics_sphere(class Brain*);
  void dump(class Brain*);

  void dump_txt(Brain*, vector<string>);
  void dump_mri(Brain*, vector<string>);
  void dump_mricrogl(Brain*, vector<string>);

  void sort_tag(Brain*, double*, int);
  tagint find_tag(Brain*, double, double, double);
  //nifti_image *nifti_image_setup(Brain*, vector<string>, int);
  nifti_image *nifti_image_setup(Brain*, vector<string>, const int[], int);

  fstream fb;
  FILE* fw;
  string rname, sname;
  bool do_dump, do_statistics;
  int revery,severy;

  vector<vector<string>> dump_arg;

};

}

#endif
