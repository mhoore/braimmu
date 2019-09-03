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

#define PI 3.141592653589793

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

#endif
