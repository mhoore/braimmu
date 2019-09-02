#ifndef SCENARIO_GEOMETRY_H
#define SCENARIO_GEOMETRY_H

#include "virtualbrain.h"

using namespace std;

class ScenarioGeometry : public VirtualBrain {
 public:
  ScenarioGeometry(int, char**, int, int);
  ~ScenarioGeometry();

  void allocations();
  void integrate(int);
  void update();
  void derivatives();
  int set_property(string,string);
  int find_agent(string);
  void set_parameters();

  ns_geometry::properties prop;

};

#endif
