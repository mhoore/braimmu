#include "ScenarioConnectomeStrategyCUDA.h"
#include "ScenarioConnectomeStrategyCUDANewton.h"

#include "scenario_connectome.h"

#include "cudaError.h"

#include <cuda.h>

static __global__ void derivativeKernelNewton(const double* agent, double* deriv, const int* type,
                                        const ScenarioConnectomeStrategyCUDANewton::array_properties arr_prop,
                                        int nall, double dt, int step, int parity);

/* ----------------------------------------------------------------------*/
void ScenarioConnectomeStrategyCUDANewton::derivatives() {

  static constexpr int BLOCK_DIM = 128;
  // set derivatives of all voxels to zero
  CUDA_SAFE_CALL( cudaMemsetAsync(deriv, 0, ScenarioConnectomeAgents::num_agents*sizeof(double)*m_this->nall) );

  //const dim3 blocks(m_this->nvl[0]/BLOCK_DIM + (m_this->nvl[0]%BLOCK_DIM>0), m_this->nvl[1], m_this->nvl[2]);
  //derivativeKernel<<<blocks, BLOCK_DIM>>>(agent, deriv, type, arr_prop, m_this->nall,m_this->dt, m_this->step);
  derivativeKernelNewton<<<m_this->nall/BLOCK_DIM + (m_this->nall%BLOCK_DIM>0), BLOCK_DIM>>>(agent, deriv, type, arr_prop, m_this->nall,m_this->dt, m_this->step, 0);
  derivativeKernelNewton<<<m_this->nall/BLOCK_DIM + (m_this->nall%BLOCK_DIM>0), BLOCK_DIM>>>(agent, deriv, type, arr_prop, m_this->nall,m_this->dt, m_this->step, 1);
}

static __device__ int findParity(Coord coord)
{
  return (int) ( (coord.z ^ coord.y ^ (coord.x % 2) ) & 1 );
  //printf("%i %i %i %i \n", coord.z, coord.y, coord.x, par);
  //return par;
}

static __global__ void derivativeKernelNewton(const double* agent, double* deriv, const int* type,
                                        const ScenarioConnectomeStrategyCUDANewton::array_properties arr_prop,
                                        int nall, double dt, int step, int parity)
{
  const int i = threadIdx.x + blockDim.x*blockIdx.x;
 
  if(i < nall) {

    Coord coord = find_coord(i);

    if (!(parity == findParity(coord))) return;

    if (type[i] & tissue(EMP)) return;

    // sAb, fAb, and tau efflux from CSF
    if (type[i] & tissue(CSF)) {
      deriv[sAb * nall + i] -= prop.es * agent[sAb * nall + i];
      deriv[fAb * nall + i] -= prop.es * agent[fAb * nall + i];
      deriv[phr * nall + i] -= prop.ephi * agent[phr * nall + i];
    }

    // in parenchyma (WM and GM)
    else {
      double dum = prop.kp * agent[sAb * nall + i] * agent[fAb * nall + i]
                     + prop.kn * agent[sAb * nall + i] * agent[sAb * nall + i];

      // sAb
      deriv[sAb * nall + i] += agent[neu * nall + i] * agent[cir * nall + i]
                            - dum
                            - prop.ds * agent[mic * nall + i] * agent[sAb * nall + i];
      // fAb
      deriv[fAb * nall + i] += dum
                            - prop.df * agent[mic * nall + i] * agent[fAb * nall + i];

      dum = prop.ktau * agent[phr * nall + i];

      // tau protein phosphorylation due to fAb and neu
      deriv[phr * nall + i] += prop.kphi * agent[fAb * nall + i] * agent[neu * nall + i]
                            - dum;

      // tau tangle formation from phosphorylated tau
      deriv[tau * nall + i] += dum;

      // neuronal death due to tau aggregation
      deriv[neu * nall + i] -= prop.dnt * agent[tau * nall + i] * agent[neu * nall + i];

      // astrogliosis
      dum = agent[fAb * nall + i] * agent[mic * nall + i];
      deriv[ast * nall + i] = prop.ka * (dum / (dum + prop.Ha) - agent[ast * nall + i]);

    }
  }
}

