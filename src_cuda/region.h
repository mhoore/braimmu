#ifndef REGION_H
#define REGION_H

#include "pointers.h"

using namespace std;

class Region {
 public:
  Region();
  ~Region();

  void apply_regions(class VirtualBrain*);
  int sphere(class VirtualBrain*, vector<string>);
  int block(class VirtualBrain*, vector<string>);

  vector<vector<string>> reg_arg;

};

#endif
