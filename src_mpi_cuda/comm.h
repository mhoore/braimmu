#ifndef COMM_H
#define COMM_H

#include <vector>
#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#include "pointers.h"
#include "memory.h"
#include "brain.h"

using namespace std;
namespace brain_NS {

class Comm : public Memory {
 public:
  Comm();
  ~Comm();

  void partition(class Brain*);
  void balance(class Brain*);
  void comm_init(class Brain*);

  void forward_comm(class Brain*);
  void forward_pack(class Brain*, int);
  void forward_unpack(class Brain*);

  void reverse_comm(class Brain*);
  void reverse_pack(class Brain*, int);
  void reverse_unpack(class Brain*);

  void allocations(class Brain*);

  int comm_side[6],buf_size[6];
  double *send_buf,*recv_buf;
  int comm_size, max_buf_size;

  int b_itr; // total number of balance iterations
  string b_dim; // balance dimension: x,y, or z

 private:
  int find_me(class Brain*, int, int, int);
  int find_id(class Brain*, int, int, int);
};

}

#endif
