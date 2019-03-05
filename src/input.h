#ifndef INPUT_H
#define INPUT_H

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
#include "memory.h"

#include "nifti1.h"
#include "nifti1_io.h"
#include "znzlib.h"

using namespace std;
namespace brain_NS {

class Input : public Memory {
 public:
  Input();
  ~Input();

  void file(const char*, class Brain*);
  void parse(string);
  void execute_command(class Brain*);

  int read_parameters(class Brain*);
  int read_mri(class Brain*);
  void read_region(class Brain*);
  int read_restart(class Brain*);
  int read_statistics(class Brain*);
  int read_dump(class Brain*);
  int find_agent(string);

  ifstream fr;

  int narg;
  string command;

  vector<string> *arg;

};

}

#endif
