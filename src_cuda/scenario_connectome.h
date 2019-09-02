#ifndef SCENARIO_CONNECTOME_H
#define SCENARIO_CONNECTOME_H

#include "virtualbrain.h"

using namespace std;

class ScenarioConnectome : public VirtualBrain {
 public:
  ScenarioConnectome(int, char**, int, int);
  ~ScenarioConnectome();

  void allocations();
  void integrate(int);
  void update();
  void derivatives();
  int set_property(string,string);
  int find_agent(string);
  void set_parameters();

  ns_connectome::properties prop;

};

#endif
