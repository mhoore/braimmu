#ifndef POINTERS_H
#define POINTERS_H

#include <vector>
#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#include <inttypes.h>

using namespace std;

namespace brain_NS {

/// in case number of voxels is more than 2 bilions
typedef int64_t tagint;
#ifndef PRId64
#define PRId64 "ld"
#endif
#define TAGINT_FORMAT "%" PRId64
/// otherwise
//#define TAGINT_FORMAT "%d"
//typedef int32_t tagint;

const string flog = "log.braimmu";

enum{XLO,XHI,YLO,YHI,ZLO,ZHI};

//// model
/// mic, ast, neu; number of microglia, atrocytes, and neurons at each voxel
/// ilb1, il6, tnf; concentration of cytokines at each voxel
/// sAb, fAb; concentration of soluble and fibrillar Amyloid beta at each voxel
enum{mic = 0,neu,sAb,fAb,ast,phr,tau,cir, num_agents}; // ast,ilb1,il6,tnf
const string ag_str[num_agents] = {"mic","neu","sAb","fAb","ast","phr","tau","cir"}; // "ilb1","il6","tnf"

/// voxel types
const int EMP_type = -1; // Empty space type
const int CSF_type = 0; // Cerebrospinal fluid type
const int WM_type = 1; // White matter parenchyma type
const int GM_type = 2; // Gray matter parenchyma type

const int ndim = 3; // number of dimensions (3D)
}

#endif
