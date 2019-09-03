#ifndef COMM_H
#define COMM_H

#include "pointers.h"

using namespace std;

class Comm {
 public:
  Comm(class VirtualBrain*);
  ~Comm();

  void partition(class VirtualBrain*);
  void balance(class VirtualBrain*);
  void comm_init(class VirtualBrain*);

  void forward_comm(class VirtualBrain*);
  void forward_pack(class VirtualBrain*, int);
  void forward_unpack(class VirtualBrain*, int);

  void reverse_comm(class VirtualBrain*);
  void reverse_pack(class VirtualBrain*, int);
  void reverse_unpack(class VirtualBrain*, int);

  void allocations(class VirtualBrain*);

  int comm_side[6];
  //double *send_buf,*recv_buf;
  vector<double> send_buf,recv_buf;
  int comm_size, max_buf_size;

  int b_itr; // total number of balance iterations
  string b_dim; // balance dimension: x,y, or z

 private:
  int find_me(class VirtualBrain*, int, int, int);

};

#endif
