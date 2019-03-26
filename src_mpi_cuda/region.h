#ifndef REGION_H
#define REGION_H

#include <vector>
#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#include "pointers.h"
#include "brain.h"

namespace brain_NS {

class Region {
 public:
  Region();
  ~Region();

  void apply_regions(class Brain*);
  int sphere(class Brain*, vector<string>);
  int block(class Brain*, vector<string>);

  vector<vector<string>> reg_arg;

};

}

#endif
