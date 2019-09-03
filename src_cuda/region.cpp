#include "virtualbrain.h"
#include "region.h"

using namespace std;

/* ---------------------------------------------------------------------- */
Region::Region() {
}

/* ----------------------------------------------------------------------*/
Region::~Region() {
  for (int i=0; i<reg_arg.size(); i++) {
    reg_arg[i].clear();
  }

  reg_arg.clear();

}

/* ----------------------------------------------------------------------
 * Apply all the regions, respectively in the order of input
 * ----------------------------------------------------------------------*/
void Region::apply_regions(VirtualBrain *brn) {
  for (int i=0; i<reg_arg.size(); i++) {
    int narg = reg_arg[i].size();

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
int Region::sphere(VirtualBrain *brn, vector<string> arg) {
  int narg = arg.size();

  int nall = brn->nall;

  auto &tissue = brn->tissue;
  auto &type = brn->type;

  auto &x = brn->x;

  double center[3];
  center[0] = stof(arg[1]);
  center[1] = stof(arg[2]);
  center[2] = stof(arg[3]);

  double radius2 = stof(arg[4]);
  radius2 *= radius2;

  bool in;
  if (!arg[5].compare("in")) in = 1;
  else if (!arg[5].compare("out")) in = 0;
  else
    return 0;

  int c = 6;
  do {
    if (c+1 >= narg)
      return 0;
    if (!arg[c].compare("type")) {
      int type_one = stoi(arg[c+1]);
      for (int i=0; i<nall; i++) {
        double delx = x[0][i] - center[0];
        double dely = x[1][i] - center[1];
        double delz = x[2][i] - center[2];

        double rsq = delx*delx + dely*dely + delz*delz;

        if (in && rsq <= radius2)
          type[i] = tissue[type_one];
        else if (!in && rsq > radius2)
          type[i] = tissue[type_one];
      }
    }

    else if (brn->find_agent(arg[c]) >= 0) {
      int ag_id = brn->find_agent(arg[c]);
      double val_one = stof(arg[c+1]);
      for (int i=0; i<nall; i++) {
        double delx = x[0][i] - center[0];
        double dely = x[1][i] - center[1];
        double delz = x[2][i] - center[2];

        double rsq = delx*delx + dely*dely + delz*delz;

        if (in && rsq <= radius2)
          brn->set_agent(ag_id,i,val_one,0);
        else if (!in && rsq > radius2)
          brn->set_agent(ag_id,i,val_one,0);
      }
    }

    else return 0;

    c += 2;
  } while (c < narg);

  return 1;

}

/* ----------------------------------------------------------------------*/
int Region::block(VirtualBrain *brn, vector<string> arg) {
  int narg = arg.size();

  int nall = brn->nall;

  auto &tissue = brn->tissue;
  auto &type = brn->type;

  auto &x = brn->x;

  double blo[3],bhi[3];
  blo[0] = stof(arg[1]);
  bhi[0] = stof(arg[2]);
  blo[1] = stof(arg[3]);
  bhi[1] = stof(arg[4]);
  blo[2] = stof(arg[5]);
  bhi[2] = stof(arg[6]);

  bool in;
  if (!arg[7].compare("in")) in = 1;
  else if (!arg[7].compare("out")) in = 0;
  else return 0;

  int c = 8;
  do {
    if (c+1 >= narg)
      return 0;
    if (!arg[c].compare("type")) {
      int type_one = stoi(arg[c+1]);
      for (int i=0; i<nall; i++) {
        if (x[0][i] >= blo[0] && x[1][i] >= blo[1] && x[2][i] >= blo[2] &&
            x[0][i] <= bhi[0] && x[1][i] <= bhi[1] && x[2][i] <= bhi[2]) {
          if (in)
            type[i] = tissue[type_one];
        } else if (!in)
          type[i] = tissue[type_one];
      }
    }

    else if (brn->find_agent(arg[c]) >= 0) {
      int ag_id = brn->find_agent(arg[c]);
      double val_one = stoi(arg[c+1]);
      for (int i=0; i<nall; i++) {
        if (x[0][i] >= blo[0] && x[1][i] >= blo[1] && x[2][i] >= blo[2] &&
            x[0][i] <= bhi[0] && x[1][i] <= bhi[1] && x[2][i] <= bhi[2]) {
          if (in)
            brn->set_agent(ag_id,i,val_one,0);
        } else if (!in)
          brn->set_agent(ag_id,i,val_one,0);
      }
    }

    else return 0;

    c += 2;
  } while (c < narg);

  return 1;

}
