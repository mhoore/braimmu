#ifndef SCENARIO_GEOMETRY_H
#define SCENARIO_GEOMETRY_H

#include "virtualbrain.h"

using namespace std;
using namespace ns_connectome;

class ScenarioGeometry : public VirtualBrain {
 public:
  ScenarioGeometry(int, char**, int, int);
  ~ScenarioGeometry();

  void allocations();
  void integrate(int);
  void update();
  void derivatives();

};

#endif
