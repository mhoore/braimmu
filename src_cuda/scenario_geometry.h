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

  int get_num_agents() { return num_agents;}
  double get_agent(int ag_id, int i) { return agent[ag_id][i];}
  double get_deriv(int ag_id, int i) { return deriv[ag_id][i];}
  string get_ag_str(int ag_id) { return ag_str[ag_id];}

  void set_agent(int ag_id, int i, double val, bool flag) {
    flag ? agent[ag_id][i] += val : agent[ag_id][i] = val;
  }

  void set_deriv(int ag_id, int i, double val, bool flag) {
    flag ? deriv[ag_id][i] += val : deriv[ag_id][i] = val;
  }

  void mri_topology(nifti_image*);
  int dump_specific(vector<string> arg) { return 0;}

  void reset();

 private:
  enum {mic = 0,neu,sAb,fAb,ast,cir, num_agents};
  const string ag_str[num_agents] = {"mic","neu","sAb","fAb","ast","cir"};

  struct properties {
    double D_sAb, diff_sAb; // diffusivity of sAb
    double D_mic, diff_mic; // diffusivity of microglia
    double cs, sens_s, cf, sens_f; // microglia chemotaxis sensitivity
    double kp, kn; // rate of polymerization and nucleation
    double ds,df; // clearance rate of sAb and fAb by microglia
    double es; // rate of sAb efflux in CSF
    double Ha; // Michaelis-Menten constant for astrogliosis
    double ka; // rate of astrogliosis
    double C_cir, c_cir, tau_cir, omega_cir; // circadian rhythm parameters
    double dna; // neuronal death rate due to astogliosis
    double dnf; // neuronal death rate due to fibrillization
  } prop;

  array<vector<double>, num_agents> agent, deriv;

  uint8_t * ptr8;
  int16_t *ptr16;
  int32_t *ptr32;
  float *ptrf;
  double *ptrd;

};

#endif
