#ifndef VIRTUALBRAIN_H
#define VIRTUALBRAIN_H

#include "pointers.h"
#include "init.h"
#include "comm.h"
#include "output.h"
#include "region.h"
#include "input.h"

using namespace std;

class VirtualBrain {
 public:
  virtual void allocations() = 0;
  virtual void integrate(int) = 0;
  virtual void update() = 0;
  virtual void derivatives() = 0;
  virtual int set_property(string,string) = 0;
  virtual int find_agent(string) = 0;
  virtual void set_parameters() = 0;

  virtual int get_num_agents() = 0;
  virtual double get_agent(int, int) = 0;
  virtual double get_deriv(int, int) = 0;
  virtual string get_ag_str(int) = 0;
  virtual void set_agent(int, int, double, bool) = 0;
  virtual void set_deriv(int, int, double, bool) = 0;

  virtual void mri_topology(nifti_image*) = 0;
  virtual int dump_specific(vector<string>) = 0;


  /* ----------------------------------------------------------------------
   * Find the local voxel id from local coordinates i,j,k
   * ----------------------------------------------------------------------*/
  int find_id(int i, int j, int k) {
    return i + (nvl[0] + 2) * (j + (nvl[1] + 2) * k);
  }

  int me,nproc;

  array<int, ndim> npart; // number of partitions in each dimension
  array<int, ndim> nv,nvl; // number of voxels in each dimension, global and local
  tagint nvoxel; // total number of voxels
  int nlocal, nghost, nall; // number of voxels for each core, local/ghost/all

  int step, Nrun, Nlog; // step, number of steps, log output

  double dt;  // timestep
  int nevery; // period of neural activity = nevery timesteps
  array<double, ndim> boxlo,boxhi,lbox; // box boundaries, and size
  array<double, ndim> xlo, xhi; // boundaries for each partition

  array<vector<double>, ndim> x;

  double vlen, vlen_1, vlen_2, vvol, vvol_1; // voxel size

  vector<int> tissue, type, group;
  vector<bool> is_loc; // local voxel = 1, ghost voxel = 0

  vector<tagint> tag; // tag of each voxel

  int num_conn_max; // maximum number of connections
  //tagint **conn; // connection tags for each voxel

  /// MRI image variables
  nifti_image *nim;

  bool newton_flux;

  MPI_Comm world;

  // classes
  Input *input;
  Init *init;
  Comm *comm;
  Output *output;
  Region *region;

  string scenario;

  vector<double> init_val;

};

#endif
