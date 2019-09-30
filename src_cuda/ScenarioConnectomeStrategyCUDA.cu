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

static __global__ void integrateKernel(const double* agent, double* deriv, const int* type, array_properties arr_prop, int nall);

ScenarioConnectomeStrategyCUDA::ScenarioConnectomeStrategyCUDA(ScenarioConnectome* pthis)
	: ScenarioConnectomeAbstractStrategy(pthis)
{
	CUDA_SAFE_CALL(
		cudaMalloc(&arr_prop.Dtau, 3*sizeof(double)*m_this->nall)
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
		cudaMemcpyToSymbol(&prop, &m_this->prop, sizeof(ScenarioConnectomeAgents::properties), cudaMemcpyHostToDevice)
		);
	CUDA_SAFE_CALL(
		cudaMemcpyToSymbol(&nvl, &m_this->nvl.data(), sizeof(int)*m_this->nvl.size(), cudaMemcpyHostToDevice)
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
			cudaMemCpy(arr_prop.Dtau+nall*a, m_this->arr_prop.Dtau[a].data(), nall*sizeof(double), cudaMemcpyHostToDevice)
			);
	}
	for(int a = 0; a < ScenarioConnectomeAgents::num_agents; ++a)
	{
		CUDA_SAFE_CALL(
			cudaMemCpy(agent+nall*a, m_this->agent[a].data(), nall*sizeof(double), cudaMemcpyHostToDevice)
			);
	}
	CUDA_SAFE_CALL(
		cudaMemCpy(type, m_this->type.data(), nall*sizeof(int), cudaMemcpyHostToDevice)
		);
}

void ScenarioConnectomeStrategyCUDA::pop()
{
	for(int a = 0; a < ScenarioConnectomeAgents::num_agents; ++a)
	{
		CUDA_SAFE_CALL(
			cudaMemCpy(m_this->agent[a].data(), agent+nall*a, nall*sizeof(double), cudaMemcpyDeviceToHost)
			);
		CUDA_SAFE_CALL(
			cudaMemCpy(m_this->deriv[a].data(), deriv+nall*a, nall*sizeof(double), cudaMemcpyDeviceToHost)
			);
	}
}

using namespace ScenarioConnectomeAgents;

/* ----------------------------------------------------------------------*/
void ScenarioConnectomeStrategyCUDA::derivatives() {

  // set derivatives of all voxels to zero
	CUDA_SAFE_CALL(
		cudaMemSet(deriv, 0, ScenarioConnectomeAgents::num_agents*sizeof(double)*m_this->nall)
		);

	static constexpr int BLOCK_DIM = 128;
	const dim3 blocks(nvl[0]/BLOCK_DIM, nvl[1], nvl[2]);
	integrateKernel<<< BLOCK_DIM, blocks>>>(agent, deriv, type, arr_prop); 
}

static __device__ int find_id(int i, int j, int k)
{
	return i + (nvl[0] + 2) * (j + (nvl[1] + 2) * k);
}

static __global__ void integrateKernel(const double* agent, double* deriv, const int* type, array_properties arr_prop)
{
	const int ii = threadIdx.x + blockDim.x*blockIdx.x +1;
	const int jj = blockIdx.y +1;
	const int kk = blockIdx.z +1;
	if(ii < nvl[0]+1)
	{
  // spatial derivatives
  /*for (int kk=1; kk<nvl[2]+1; kk++)*/
  /*  for (int jj=1; jj<nvl[1]+1; jj++)*/
  /*    for (int ii=1; ii<nvl[0]+1; ii++) {*/
        const int i = find_id(ii,jj,kk);
        if (type[i] & tissue(EMP)) continue;

        // direct function or time derivatives

        // sAb, fAb, and tau efflux from CSF
        if (type[i] & tissue(CSF)) {
          deriv[sAb][i] -= prop.es * agent[sAb][i];
          deriv[fAb][i] -= prop.es * agent[fAb][i];
          deriv[phr][i] -= prop.ephi * agent[phr][i];
        }
        // in parenchyma (WM and GM)
        else {
          double dum = prop.kp * agent[sAb][i] * agent[fAb][i]
                     + prop.kn * agent[sAb][i] * agent[sAb][i];

          // sAb
          deriv[sAb][i] += agent[neu][i] * agent[cir][i]
                            - dum
                            - prop.ds * agent[mic][i] * agent[sAb][i];
          // fAb
          deriv[fAb][i] += dum
                           - prop.df * agent[mic][i] * agent[fAb][i];

          dum = prop.ktau * agent[phr][i];

          // tau protein phosphorylation due to fAb and neu
          deriv[phr][i] += prop.kphi * agent[fAb][i] * agent[neu][i]
                         - dum;

          // tau tangle formation from phosphorylated tau
          deriv[tau][i] += dum;

          // neuronal death due to tau aggregation
          deriv[neu][i] -= prop.dnt * agent[tau][i] * agent[neu][i];

          // astrogliosis
          dum = agent[fAb][i] * agent[mic][i];
          deriv[ast][i] = prop.ka * (dum / (dum + prop.Ha) - agent[ast][i]);

          // circadian rhythm
          if (prop.c_cir > 0)
            deriv[cir][i] = - prop.C_cir * prop.c_cir * prop.omega_cir
                            * sin(prop.omega_cir * m_this->dt * m_this->step);
        }

        // spatial derivatives: fluxes
        int n_ngh, ngh[6];
        ngh[0] = m_this->find_id(ii+1,jj,kk);
        ngh[1] = m_this->find_id(ii,jj+1,kk);
        ngh[2] = m_this->find_id(ii,jj,kk+1);
        n_ngh = 3;

        /*if (!newton_flux) {*/
          ngh[3] = m_this->find_id(ii-1,jj,kk);
          ngh[4] = m_this->find_id(ii,jj-1,kk);
          ngh[5] = m_this->find_id(ii,jj,kk-1);
          n_ngh = 6;
        /*}*/

        for (int c=0; c<n_ngh; ++c) {
          int j = ngh[c];
          int d = c;
          if (c >= 3)
            d = c - 3;

          if (type[j] & tissue(EMP)) continue;

          double del_phr = agent[phr][i] - agent[phr][j];

          // diffusion of tau
          double dum = 0.5 * (arr_prop.Dtau[d][i] + arr_prop.Dtau[d][j]) * del_phr;
          deriv[phr][i] -= dum;
          /*if (newton_flux)*/
          /*  deriv[phr][j] += dum;*/

          double del_sAb = agent[sAb][i] - agent[sAb][j];

          // diffusion of sAb
          dum = prop.D_sAb * del_sAb;
          deriv[sAb][i] -= dum;
          /*if (newton_flux)*/
          /*  deriv[sAb][j] += dum;*/

          // only in parenchyma
          if (type[i] & tissue(WM) || type[i] & tissue(GM))
            if (type[j] & tissue(WM) || type[j] & tissue(GM)) {
              double del_fAb = agent[fAb][i] - agent[fAb][j];
              double del_mic = agent[mic][i] - agent[mic][j];

              // migration of microglia toward higher sAb concentrations
              dum = prop.cs * del_sAb;
              if (del_sAb > 0.0)
                dum *= agent[mic][j];
              else
                dum *= agent[mic][i];

              deriv[mic][i] += dum;
              /*if (newton_flux)*/
              /*  deriv[mic][j] -= dum;*/

              // migration of microglia toward higher fAb concentrations
              dum = prop.cf * del_fAb;
              if (del_fAb > 0.0)
                dum *= agent[mic][j];
              else
                dum *= agent[mic][i];

              deriv[mic][i] += dum;
              /*if (newton_flux)*/
              /*  deriv[mic][j] -= dum;*/

              // diffusion of microglia
              dum = prop.D_mic * del_mic;
              deriv[mic][i] -= dum;
              /*if (newton_flux)*/
              /*  deriv[mic][j] += dum;*/
            }
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
