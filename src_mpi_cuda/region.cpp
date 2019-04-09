#include <mpi.h>
#include <vector>
#include "math.h"

#include "pointers.h"
#include "brain.h"
#include "region.h"

using namespace std;
using namespace brain_NS;

/* ---------------------------------------------------------------------- */
Region::Region() {
}

/* ----------------------------------------------------------------------*/
Region::~Region() {
  int i;

  for (i=0; i<reg_arg.size(); i++) {
    reg_arg[i].clear();
  }

  reg_arg.clear();

}

/* ----------------------------------------------------------------------
 * Apply all the regions, respectively in the order of input
 * ----------------------------------------------------------------------*/
void Region::apply_regions(Brain *brn) {
  int i,j;
  int narg;

  for (i=0; i<reg_arg.size(); i++) {
    narg = reg_arg[i].size();

    if (!reg_arg[i][0].compare("sphere")) {
      if (narg < 8) {
        printf("Error: region sphere has fewer argument than required. \n");
        exit(1);
      }

      if (!sphere(brn,reg_arg[i])) {
        printf("Error: region sphere cannot be assigned. \n");
        exit(1);
      }
    }

    else if (!reg_arg[i][0].compare("block")) {
      if (narg < 10) {
        printf("Error: region block has fewer argument than required. \n");
        exit(1);
      }

      if (!block(brn,reg_arg[i])) {
        printf("Error: region block cannot be assigned. \n");
        exit(1);
      }
    }

    else {
      printf("Error: region unknown. \n");
      exit(1);
    }

  }

}

/* ----------------------------------------------------------------------*/
int Region::sphere(Brain *brn, vector<string> arg) {
  int i,c;
  double center[3],radius2,delx,dely,delz,rsq;
  bool in;

  int narg = arg.size();

  int nall = brn->nall;

  int *type = brn->type;

  double **x = brn->x;

  center[0] = stof(arg[1]);
  center[1] = stof(arg[2]);
  center[2] = stof(arg[3]);
  radius2 = stof(arg[4]);
  radius2 *= radius2;

  if (!arg[5].compare("in")) in = 1;
  else if (!arg[5].compare("out")) in = 0;
  else
    return 0;

  c = 6;
  do {
    if (c+1 >= narg)
      return 0;
    if (!arg[c].compare("type")) {
      int type_one = stoi(arg[c+1]);
      for (i=0; i<nall; i++) {
        delx = x[i][0] - center[0];
        dely = x[i][1] - center[1];
        delz = x[i][2] - center[2];

        rsq = delx*delx + dely*dely + delz*delz;

        if (in && rsq <= radius2)
          type[i] = type_one;
        else if (!in && rsq > radius2)
          type[i] = type_one;
      }
    }

    else if (brn->input->find_agent(arg[c]) >= 0) {
      int ag_id = brn->input->find_agent(arg[c]);
      double val_one = stof(arg[c+1]);
      for (i=0; i<nall; i++) {
        delx = x[i][0] - center[0];
        dely = x[i][1] - center[1];
        delz = x[i][2] - center[2];

        rsq = delx*delx + dely*dely + delz*delz;

        if (in && rsq <= radius2)
          brn->agent[ag_id][i] = val_one;
        else if (!in && rsq > radius2)
          brn->agent[ag_id][i] = val_one;
      }
    }

    else return 0;

    c += 2;
  } while (c < narg);

  return 1;

}

/* ----------------------------------------------------------------------*/
int Region::block(Brain *brn, vector<string> arg) {
  int i,c;
  double blo[3],bhi[3];
  bool in;

  int narg = arg.size();

  int nall = brn->nall;

  int *type = brn->type;

  double **x = brn->x;

  blo[0] = stof(arg[1]);
  bhi[0] = stof(arg[2]);
  blo[1] = stof(arg[3]);
  bhi[1] = stof(arg[4]);
  blo[2] = stof(arg[5]);
  bhi[2] = stof(arg[6]);

  if (!arg[7].compare("in")) in = 1;
  else if (!arg[7].compare("out")) in = 0;
  else return 0;

  c = 8;
  do {
    if (c+1 >= narg)
      return 0;
    if (!arg[c].compare("type")) {
      int type_one = stoi(arg[c+1]);
      for (i=0; i<nall; i++) {
        if (x[i][0] >= blo[0] && x[i][1] >= blo[1] && x[i][2] >= blo[2] &&
            x[i][0] <= bhi[0] && x[i][1] <= bhi[1] && x[i][2] <= bhi[2]) {
          if (in)
            type[i] = type_one;
        } else if (!in)
          type[i] = type_one;
      }
    }

    else if (brn->input->find_agent(arg[c]) >= 0) {
      int ag_id = brn->input->find_agent(arg[c]);
      double val_one = stoi(arg[c+1]);
      for (i=0; i<nall; i++) {
        if (x[i][0] >= blo[0] && x[i][1] >= blo[1] && x[i][2] >= blo[2] &&
            x[i][0] <= bhi[0] && x[i][1] <= bhi[1] && x[i][2] <= bhi[2]) {
          if (in)
            brn->agent[ag_id][i] = val_one;
        } else if (!in)
          brn->agent[ag_id][i] = val_one;
      }
    }

    else return 0;

    c += 2;
  } while (c < narg);

  return 1;

}
