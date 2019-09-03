#ifndef INPUT_H
#define INPUT_H

#include "pointers.h"

using namespace std;

class Input {
 public:
  Input();
  ~Input();

  void file(const char*, class VirtualBrain*);
  void parse(string);
  void execute_command(class VirtualBrain*);

  int read_parameters(class VirtualBrain*);
  void read_mri(class VirtualBrain*);
  void read_region(class VirtualBrain*);
  int read_restart(class VirtualBrain*);
  int read_statistics(class VirtualBrain*);
  int read_dump(class VirtualBrain*);

  ifstream fr;

  int narg;
  string command;

  vector<string> arg;

};

#endif
