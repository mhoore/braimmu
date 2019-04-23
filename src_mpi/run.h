#ifndef RUN_H
#define RUN_H

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

class Run {
 public:
  Run();
  ~Run();

  void integrate(class Brain*, int);
  void update(class Brain*);
  void derivatives(Brain *brn);

  int find_id(class Brain*, int, int, int);

};

}

#endif
