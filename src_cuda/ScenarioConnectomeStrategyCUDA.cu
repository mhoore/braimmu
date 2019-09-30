#include "ScenarioConnectomeStrategyCUDA.h"

#include "scenario_connectome.h"

#include "cudaError.h"

#include <cuda.h>

using namespace std;

__constant__ ScenarioConnectome::properties prop;
__constant__ double nvl[ndim];

static __device__ constexpr int tissue(int type)
{
	return 1<<type;
}

static __global__ void integrateKernel(const double* agent, double* deriv, const int* type, const ScenarioConnectomeStrategyCUDA::array_properties arr_prop, int nalli,
double dt, int step);

ScenarioConnectomeStrategyCUDA::ScenarioConnectomeStrategyCUDA(ScenarioConnectome* pthis)
	: ScenarioConnectomeAbstractStrategy(pthis)
{
	CUDA_SAFE_CALL(
		cudaMalloc(&arr_prop.Dtau, ndim*sizeof(double)*m_this->nall)
		);
	CUDA_SAFE_CALL(
		cudaMalloc(&agent, ScenarioConnectomeAgents::num_agents*sizeof(double)*m_this->nall)
		);
	CUDA_SAFE_CALL(
		cudaMalloc(&deriv, ScenarioConnectomeAgents::num_agents*sizeof(double)*m_this->nall)
		);
	CUDA_SAFE_CALL(
		cudaMalloc(&type, sizeof(int)*m_this->nall)
		);

	CUDA_SAFE_CALL(
		cudaMemcpyToSymbol(&prop, &m_this->prop, sizeof(ScenarioConnectome::properties), cudaMemcpyHostToDevice)
		);
	CUDA_SAFE_CALL(
		cudaMemcpyToSymbol(&nvl, m_this->nvl.data(), sizeof(int)*m_this->nvl.size(), cudaMemcpyHostToDevice)
		);
}

ScenarioConnectomeStrategyCUDA::~ScenarioConnectomeStrategyCUDA()
{
	CUDA_SAFE_CALL(
		cudaFree(arr_prop.Dtau)
		);
	CUDA_SAFE_CALL(
		cudaFree(agent)
		);
	CUDA_SAFE_CALL(
		cudaFree(deriv)
		);
	CUDA_SAFE_CALL(
		cudaFree(type)
		);
}

void ScenarioConnectomeStrategyCUDA::push()
{
	for(int a = 0; a < ndim; ++a)
	{
		CUDA_SAFE_CALL(
			cudaMemcpy(arr_prop.Dtau+m_this->nall*a, m_this->arr_prop.Dtau[a].data(), m_this->nall*sizeof(double), cudaMemcpyHostToDevice)
			);
	}
	for(int a = 0; a < ScenarioConnectomeAgents::num_agents; ++a)
	{
		CUDA_SAFE_CALL(
			cudaMemcpy(agent+m_this->nall*a, m_this->agent[a].data(), m_this->nall*sizeof(double), cudaMemcpyHostToDevice)
			);
	}
	CUDA_SAFE_CALL(
		cudaMemcpy(type, m_this->type.data(), m_this->nall*sizeof(int), cudaMemcpyHostToDevice)
		);
}

void ScenarioConnectomeStrategyCUDA::pop()
{
	for(int a = 0; a < ScenarioConnectomeAgents::num_agents; ++a)
	{
		CUDA_SAFE_CALL(
			cudaMemcpy(m_this->agent[a].data(), agent + m_this->nall*a, m_this->nall*sizeof(double), cudaMemcpyDeviceToHost)
			);
		CUDA_SAFE_CALL(
			cudaMemcpy(m_this->deriv[a].data(), deriv+m_this->nall*a, m_this->nall*sizeof(double), cudaMemcpyDeviceToHost)
			);
	}
}

using namespace ScenarioConnectomeAgents;

/* ----------------------------------------------------------------------*/
void ScenarioConnectomeStrategyCUDA::derivatives() {

  // set derivatives of all voxels to zero
	CUDA_SAFE_CALL(
		cudaMemset(deriv, 0, ScenarioConnectomeAgents::num_agents*sizeof(double)*m_this->nall)
		);

	static constexpr int BLOCK_DIM = 128;
	const dim3 blocks(m_this->nvl[0]/BLOCK_DIM, m_this->nvl[1], m_this->nvl[2]);
	integrateKernel<<< BLOCK_DIM, blocks>>>(agent, deriv, type, arr_prop, m_this->nall,m_this->dt, m_this->step); 
}

static __device__ int find_id(int i, int j, int k)
{
	return i + (nvl[0] + 2) * (j + (nvl[1] + 2) * k);
}

static __global__ void integrateKernel(const double* agent, double* deriv, const int* type, const ScenarioConnectomeStrategyCUDA::array_properties arr_prop, int nall,
double dt, int step)
{
	const int ii = threadIdx.x + blockDim.x*blockIdx.x +1;
	const int jj = blockIdx.y +1;
	const int kk = blockIdx.z +1;
	if(ii < nvl[0]+1)
	{
        const int i = find_id(ii,jj,kk);
        if (type[i] & tissue(EMP)) return;

        // direct function or time derivatives

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

        // spatial derivatives: fluxes
        int n_ngh, ngh[6];
        ngh[0] = find_id(ii+1,jj,kk);
        ngh[1] = find_id(ii,jj+1,kk);
        ngh[2] = find_id(ii,jj,kk+1);
        n_ngh = 3;

        /*if (!newton_flux) {*/
          ngh[3] = find_id(ii-1,jj,kk);
          ngh[4] = find_id(ii,jj-1,kk);
          ngh[5] = find_id(ii,jj,kk-1);
          n_ngh = 6;
        /*}*/

        for (int c=0; c<n_ngh; ++c) {
          int j = ngh[c];
          int d = c;
          if (c >= 3)
            d = c - 3;

          if (type[j] & tissue(EMP)) continue;

          double del_phr = agent[phr * nall + i] - agent[phr * nall + j];

          // diffusion of tau
          double dum = 0.5 * (arr_prop.Dtau[ ndim * i + d] + arr_prop.Dtau[ndim * j + d]) * del_phr;
          deriv[phr * nall + i] -= dum;
          /*if (newton_flux)*/
          /*  deriv[phr][j] += dum;*/

          double del_sAb = agent[sAb * nall + i] - agent[sAb * nall + j];

          // diffusion of sAb
          dum = prop.D_sAb * del_sAb;
          deriv[sAb * nall + i] -= dum;
          /*if (newton_flux)*/
          /*  deriv[sAb][j] += dum;*/

          // only in parenchyma
          if (type[i] & tissue(WM) || type[i] & tissue(GM))
            if (type[j] & tissue(WM) || type[j] & tissue(GM)) {
              double del_fAb = agent[fAb * nall + i] - agent[fAb * nall + j];
              double del_mic = agent[mic * nall + i] - agent[mic * nall + j];

              // migration of microglia toward higher sAb concentrations
              dum = prop.cs * del_sAb;
              if (del_sAb > 0.0)
                dum *= agent[mic * nall + j];
              else
                dum *= agent[mic * nall + i];

              deriv[mic * nall + i] += dum;
              /*if (newton_flux)*/
              /*  deriv[mic][j] -= dum;*/

              // migration of microglia toward higher fAb concentrations
              dum = prop.cf * del_fAb;
              if (del_fAb > 0.0)
                dum *= agent[mic * nall + j];
              else
                dum *= agent[mic * nall + i];

              deriv[mic * nall + i] += dum;
              /*if (newton_flux)*/
              /*  deriv[mic][j] -= dum;*/

              // diffusion of microglia
              dum = prop.D_mic * del_mic;
              deriv[mic * nall + i] -= dum;
              /*if (newton_flux)*/
              /*  deriv[mic][j] += dum;*/
      }
    }
	}
}

/* ----------------------------------------------------------------------*/
void ScenarioConnectomeStrategyCUDA::update() {

	pop();
	using namespace ScenarioConnectomeAgents;

  // update local voxels
  for (int kk=1; kk<m_this->nvl[2]+1; kk++)
    for (int jj=1; jj<m_this->nvl[1]+1; jj++)
      for (int ii=1; ii<m_this->nvl[0]+1; ii++) {
        int i = m_this->find_id(ii,jj,kk);
        if (m_this->type[i] & m_this->tissue[EMP]) continue;

        // time integration (Euler's scheme)
        for (int ag_id=0; ag_id<num_agents; ag_id++)
          m_this->agent[ag_id][i] += m_this->deriv[ag_id][i] * m_this->dt;
      }

	push();
}
