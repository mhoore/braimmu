#ifndef SCENARIO_CONNECTOME_H
#define SCENARIO_CONNECTOME_H

#include "virtualbrain.h"

using namespace std;
using namespace ns_connectome;

class ScenarioConnectome : public VirtualBrain {
 public:
  ScenarioConnectome(int, char**, int, int);
  ~ScenarioConnectome();

  void allocations();
  void integrate(int);
  void update();
  void derivatives();

};

#endif
