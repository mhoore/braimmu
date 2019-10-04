#include "ScenarioConnectomeStrategyCUDA.h"
#include "ScenarioConnectomeStrategyCUDANewton.h"

#include "scenario_connectome.h"

#include "cudaError.h"

#include <cuda.h>

static __global__ void derivativeKernelNewton(const double* agent, double* deriv, const int* type,
                                        const ScenarioConnectomeStrategyCUDA::array_properties arr_prop,
                                        int nall, double dt, int step, int parity);

/* ----------------------------------------------------------------------*/
void ScenarioConnectomeStrategyCUDANewton::derivatives() {

  static constexpr int BLOCK_DIM = 128;
  // set derivatives of all voxels to zero
  CUDA_SAFE_CALL( cudaMemsetAsync(deriv, 0, ScenarioConnectomeAgents::num_agents*sizeof(double)*m_this->nall) );

  const dim3 blocks( (m_this->nvl[0] + 1)/2/BLOCK_DIM + (( (m_this->nvl[0] + 1)/2)%BLOCK_DIM>0), m_this->nvl[1], m_this->nvl[2]);
  derivativeKernelNewton<<<blocks, BLOCK_DIM>>>(agent, deriv, type, arr_prop, m_this->nall,m_this->dt, m_this->step, 0);
  derivativeKernelNewton<<<blocks, BLOCK_DIM>>>(agent, deriv, type, arr_prop, m_this->nall,m_this->dt, m_this->step, 1);
  //derivativeKernelNewton<<<m_this->nall/BLOCK_DIM + (m_this->nall%BLOCK_DIM>0), BLOCK_DIM>>>(agent, deriv, type, arr_prop, m_this->nall,m_this->dt, m_this->step, 0);
  //derivativeKernelNewton<<<m_this->nall/BLOCK_DIM + (m_this->nall%BLOCK_DIM>0), BLOCK_DIM>>>(agent, deriv, type, arr_prop, m_this->nall,m_this->dt, m_this->step, 1);
}

static __device__ int findParity(Coord coord)
{
  return (int) ( (coord.z ^ coord.y ^ (coord.x % 2) ) & 1 );
}

static __global__ void derivativeKernelNewton(const double* agent, double* deriv, const int* type,
                                        const ScenarioConnectomeStrategyCUDA::array_properties arr_prop,
                                        int nall, double dt, int step, int parity)
{
  Coord coord;

  coord.x = threadIdx.x + blockDim.x*blockIdx.x;
  coord.y = blockIdx.y +1;
  coord.z = blockIdx.z +1;

  // choose the right parity
  coord.x = (( (coord.z ^ coord.y) ^ parity) & 1) + 2 * coord.x;
  
  if(coord.x < nvl[0]+1 && coord.x > 0) {
    const int i = find_id(coord.x,coord.y,coord.z);

    //printf("id=%li , coords = %i %i %i parity = %i \n", i, coord.x, coord.y, coord.z, findParity(coord));

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

      // circadian rhythm
      if (prop.c_cir > 0)
        deriv[cir * nall + i] = - prop.C_cir * prop.c_cir * prop.omega_cir
                            * sin(prop.omega_cir * dt * step);
    }

    for (int d=0; d < 3; d+=1) {
		  const int j = find_id(coord.x + (d==0), coord.y + (d==1), coord.z + (d==2));

		  if (type[j] & tissue(EMP)) continue;

		  double del_phr = agent[phr * nall + i] - agent[phr * nall + j];

		  // diffusion of tau
      double dum = 0.5 * (arr_prop.Dtau[ nall * d + i] + arr_prop.Dtau[nall * d + j]) * del_phr;
      deriv[phr * nall + i] -= dum;
      deriv[phr * nall + j] += dum;

		  double del_sAb = agent[sAb * nall + i] - agent[sAb * nall + j];

		  // diffusion of sAb
		  dum = prop.D_sAb * del_sAb;
      deriv[sAb * nall + i] -= dum;
      deriv[sAb * nall + j] += dum;

      //if (j == 30)
      //printf("i=%li , j= %li, sAb = %g \n", i, j, deriv[sAb * nall + j]);

		  // only in parenchyma
		  if (type[i] & tissue(WM) || type[i] & tissue(GM))
			if (type[j] & tissue(WM) || type[j] & tissue(GM)) {
			  double del_fAb = agent[fAb * nall + i] - agent[fAb * nall + j];
			  double del_mic = agent[mic * nall + i] - agent[mic * nall + j];

			  // migration of microglia toward higher sAb concentrations
			  dum = prop.cs * del_sAb * agent[mic * nall + ((del_sAb > 0.0) ? j : i)];
        deriv[mic * nall + i] += dum;
        deriv[mic * nall + j] -= dum;

			  // migration of microglia toward higher fAb concentrations
        dum = prop.cf * del_fAb * agent[mic * nall + ((del_fAb > 0.0) ? j : i)];
        deriv[mic * nall + i] += dum;
        deriv[mic * nall + j] -= dum;

			  // diffusion of microglia
			  dum = prop.D_mic * del_mic;
        deriv[mic * nall + i] -= dum;
        deriv[mic * nall + j] += dum;
		  }
    }
  }
}

/* ----------------------------------------------------------------------*/
void ScenarioConnectomeStrategyCUDANewton::update() {

  using namespace ScenarioConnectomeAgents;

	static constexpr int BLOCK_DIM = 128;
	//const dim3 blocks(m_this->nvl[0]/BLOCK_DIM + (m_this->nvl[0]%BLOCK_DIM>0), m_this->nvl[1], m_this->nvl[2]);
	updateKernel<<<m_this->nall/BLOCK_DIM + (m_this->nall%BLOCK_DIM>0), BLOCK_DIM>>>(agent, deriv, type, arr_prop,m_this->dt, m_this->nall);

}

static __global__ void updateKernel(double* agent, const double* deriv, const int* type, const ScenarioConnectomeStrategyCUDA::array_properties arr_prop, double dt, int nall)
{
	const int i = threadIdx.x + blockDim.x*blockIdx.x;
	//const int jj = blockIdx.y +1;
	//const int kk = blockIdx.z +1;
  Coord coord = find_coord(i);
  if (coord.x > 0 && coord.x <= nvl[0]
   && coord.y > 0 && coord.y <= nvl[1]
   && coord.z > 0 && coord.z <= nvl[2]) {
    if (type[i] & tissue(EMP)) return;
    // time integration (Euler's scheme)
    for (int ag_id=0; ag_id<ScenarioConnectomeAgents::num_agents; ag_id++)
      agent[ag_id * nall + i] += deriv[ag_id * nall + i] * dt;
  }
}
