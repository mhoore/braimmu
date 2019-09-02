#ifndef POINTERS_H
#define POINTERS_H

#include <mpi.h>
#include <cmath>
#include <vector>

#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#include "nifti1.h"
#include "nifti1_io.h"
#include "znzlib.h"

#include <stdint.h>

#include <inttypes.h>

using namespace std;

/// in case number of voxels is more than 2 bilions
typedef int64_t tagint;
#ifndef PRId64
#define PRId64 "ld"
#endif
#define TAGINT_FORMAT "%" PRId64
/// otherwise
//#define TAGINT_FORMAT "%d"
//typedef int32_t tagint;

// union data struct for packing 32-bit and 64-bit ints into double bufs
// this avoids aliasing issues by having 2 pointers (double,int)
union ubuf {
  double d;
  int64_t i;
  ubuf(double arg) : d(arg) {}
  ubuf(int64_t arg) : i(arg) {}
  ubuf(int arg) : i(arg) {}
};

const string flog = "log.braimmu";

enum{XLO,XHI,YLO,YHI,ZLO,ZHI};

/// voxel tissue types
/* EMP : Empty space
 * CSF : Cerebrospinal fluid
 * WM : White matter parenchyma
 * GM : Grey matter parenchyma */
enum{EMP = 0,CSF,WM,GM, num_types};

const int ndim = 3; // number of dimensions

//// model
enum {mic = 0,neu,sAb,fAb,ast,phr,tau,cir, num_agents};
const string ag_str[num_agents] = {"mic","neu","sAb","fAb","ast","phr","tau","cir"};

namespace ns_connectome {
  struct properties {
    array<vector<double>, ndim> Dtau; // diffusion tensor for tau protein
    double Dtau_max, diff_tau; // maximum diffusion of tau protein
    double dnt; // neuronal death rate due to tau accumulation
    double D_sAb, diff_sAb; // diffusivity of sAb
    double D_mic, diff_mic; // diffusivity of microglia
    double cs, sens_s, cf, sens_f; // microglia chemotaxis sensitivity
    double kp, kn; // rate of polymerization and nucleation
    double ds,df; // clearance rate of sAb and fAb by microglia
    double es; // rate of sAb efflux in CSF
    double Ha; // Michaelis-Menten constant for astrogliosis
    double ka; // rate of astrogliosis
    double C_cir, c_cir, tau_cir, omega_cir; // circadian rhythm parameters
    double ktau; // rate of tau tangle formation from phosphorylated tau
    double kphi; // phosphorylation rate of tau proteins due to F and N
    double ephi; // rate of phosphorylated-tau efflux in CSF

    double dna; // neuronal death rate due to astogliosis
    double dnf; // neuronal death rate due to fibrillization
  };
}

namespace ns_geometry {
/*  //// model
  enum {mic = 0,neu,sAb,fAb,ast,cir, num_agents};
  const string ag_str[num_agents] = {"mic","neu","sAb","fAb","ast","cir"};
*/
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
  };

}

#endif
