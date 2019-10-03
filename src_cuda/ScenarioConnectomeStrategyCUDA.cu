#include "ScenarioConnectomeStrategyCUDA.h"

#include "scenario_connectome.h"

#include "cudaError.h"

#include <cuda.h>

using namespace std;

__constant__ ScenarioConnectome::properties prop;
__constant__ int nvl[ndim];
__constant__ ScenarioConnectomeStrategyCUDA::AllocPitch pitch;

static __device__ constexpr int tissue(int type)
{
	return 1<<type;
}

static __global__ void derivativeKernel(const double* agent, double* deriv, const int* type, const ScenarioConnectomeStrategyCUDA::array_properties arr_prop, int nall,
double dt, int step);

static __global__ void updateKernel(double* agent, const double* deriv, const int* type, const ScenarioConnectomeStrategyCUDA::array_properties arr_prop, double dt, int nall);

size_t ScenarioConnectomeStrategyCUDA::nvx()
{
	return m_this->nvl[0]+2;
}
size_t ScenarioConnectomeStrategyCUDA::nvyz()
{
	return (m_this->nvl[1]+2)*(m_this->nvl[2]+2);
}

ScenarioConnectomeStrategyCUDA::ScenarioConnectomeStrategyCUDA(ScenarioConnectome* pthis)
	: ScenarioConnectomeAbstractStrategy(pthis)
{
	const size_t width = this->nvx();
	const size_t height = this->nvyz();
	CUDA_SAFE_CALL(
		cudaMallocPitch(&arr_prop.Dtau, &m_allocPitch.pDouble, width*sizeof(double), height*ndim)
		);
	std::cout << "pitch Dtau: " << m_allocPitch.pDouble << std::endl;
	CUDA_SAFE_CALL(
		cudaMallocPitch(&agent, &m_allocPitch.pDouble, width*sizeof(double), height*ScenarioConnectomeAgents::num_agents)
		);
	std::cout << "pitch agent: " << m_allocPitch.pDouble << std::endl;
	CUDA_SAFE_CALL(
		cudaMallocPitch(&deriv, &m_allocPitch.pDouble, width*sizeof(double), height*ScenarioConnectomeAgents::num_agents)
		);
	std::cout << "pitch deriv: " << m_allocPitch.pDouble << std::endl;
	CUDA_SAFE_CALL(
		cudaMallocPitch(&type, &m_allocPitch.pInt, width*sizeof(int), height)
		);
	std::cout << "pitch type: " << m_allocPitch.pInt << std::endl;

	m_allocPitch.pDouble /= sizeof(double);
	m_allocPitch.pInt /= sizeof(int);

	CUDA_SAFE_CALL(
		cudaMemcpyToSymbol(prop, &m_this->prop, sizeof(ScenarioConnectome::properties))
		);
	CUDA_SAFE_CALL(
		cudaMemcpyToSymbol(nvl, (void*)m_this->nvl.data(), sizeof(int)*m_this->nvl.size())
		);
	CUDA_SAFE_CALL(
		cudaMemcpyToSymbol(pitch, &m_allocPitch, sizeof(AllocPitch))
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
	const size_t width = this->nvx();
	const size_t height = this->nvyz();
	for(int a = 0; a < ndim; ++a)
	{
		CUDA_SAFE_CALL(
			/*cudaMemcpy(arr_prop.Dtau+m_this->nall*a, m_this->arr_prop.Dtau[a].data(), m_this->nall*sizeof(double), cudaMemcpyHostToDevice)*/
			cudaMemcpy2D(arr_prop.Dtau+height*m_allocPitch.pDouble*a, m_allocPitch.pDouble*sizeof(double)
				, m_this->arr_prop.Dtau[a].data(), width*sizeof(double), width*sizeof(double), height
				, cudaMemcpyHostToDevice)
			);
	}
	for(int a = 0; a < ScenarioConnectomeAgents::num_agents; ++a)
	{
		CUDA_SAFE_CALL(
			/*cudaMemcpy(agent+m_this->nall*a, m_this->agent[a].data(), m_this->nall*sizeof(double), cudaMemcpyHostToDevice)*/
			cudaMemcpy2D(agent+height*m_allocPitch.pDouble*a, m_allocPitch.pDouble*sizeof(double)
				, m_this->agent[a].data(), width*sizeof(double), width*sizeof(double), height
				, cudaMemcpyHostToDevice)
			);
	}
	CUDA_SAFE_CALL(
		cudaMemsetAsync(type, 0, sizeof(int)*m_allocPitch.pInt*height)
		);
	CUDA_SAFE_CALL(
		/*cudaMemcpy(type, m_this->type.data(), m_this->nall*sizeof(int), cudaMemcpyHostToDevice)*/
		cudaMemcpy2D(type, m_allocPitch.pInt*sizeof(int)
			, m_this->type.data(), width*sizeof(int), width*sizeof(int), height
			, cudaMemcpyHostToDevice)
		);
}

void ScenarioConnectomeStrategyCUDA::pop()
{
	int a = m_this->find_id(26,94,75);
	printf("before pop: agent[tau] %g %d\n", (float)m_this->agent[tau][a], 23);
	printf("before pop: agent[tau] %g %d\n", m_this->agent[tau][a], 23);
	std::cout << "before pop: agent[tau][25] " << m_this->agent[tau][a] << " " << a << std::endl;
	const size_t width = this->nvx();
	const size_t height = this->nvyz();
	for(int a = 0; a < ScenarioConnectomeAgents::num_agents; ++a)
	{
		CUDA_SAFE_CALL(
			/*cudaMemcpy(m_this->agent[a].data(), agent + m_this->nall*a, m_this->nall*sizeof(double), cudaMemcpyDeviceToHost)*/
			cudaMemcpy2D(m_this->agent[a].data(), width*sizeof(double)
				, agent+height*m_allocPitch.pDouble*a, m_allocPitch.pDouble*sizeof(double)
				, width*sizeof(double), height
				, cudaMemcpyDeviceToHost)
			);
		CUDA_SAFE_CALL(
			/*cudaMemcpy(m_this->deriv[a].data(), deriv+m_this->nall*a, m_this->nall*sizeof(double), cudaMemcpyDeviceToHost)*/
			cudaMemcpy2D(m_this->deriv[a].data(), width*sizeof(double)
				, deriv+height*m_allocPitch.pDouble*a, m_allocPitch.pDouble*sizeof(double)
				, width*sizeof(double), height
				, cudaMemcpyDeviceToHost)
			);
	}

	std::cout << "after pop: agent[tau][25] " << m_this->agent[tau][a] << " " << a << std::endl;
}

using namespace ScenarioConnectomeAgents;

/*static __global__ void zeroKernel(double* data, size_t n)
{
	int* tmp = (int*)data;
	const int id = threadIdx.x + blockIdx.x*blockDim.x;
	tmp[id] = 0;
	tmp[id+n] = 0;
}*/

/* ----------------------------------------------------------------------*/
void ScenarioConnectomeStrategyCUDA::derivatives() {

	static constexpr int BLOCK_DIM = 128;
  // set derivatives of all voxels to zero
#if 0
	{
		const size_t n = ScenarioConnectomeAgents::num_agents*m_this->nall;
		zeroKernel<<<n/BLOCK_DIM + (n%BLOCK_DIM>0), BLOCK_DIM>>>(deriv, n);
	}
#endif
	CUDA_SAFE_CALL(
		cudaMemsetAsync(deriv, 0, ScenarioConnectomeAgents::num_agents*sizeof(double)*m_allocPitch.pDouble*nvyz())
	);

	const dim3 blocks(m_this->nvl[0]/BLOCK_DIM + (m_this->nvl[0]%BLOCK_DIM>0), m_this->nvl[1], m_this->nvl[2]);
	derivativeKernel<<<blocks, BLOCK_DIM>>>(agent, deriv, type, arr_prop, m_this->nall,m_this->dt, m_this->step);
}

static __device__ int find_id(int i, int j, int k)
{
	return i + (pitch.pDouble) * (j + (nvl[1] + 2) * k);
}

static __global__ void derivativeKernel(const double* agent, double* deriv, const int* type, const ScenarioConnectomeStrategyCUDA::array_properties arr_prop, int pnall,
double dt, int step)
{
	const int ii = threadIdx.x + blockDim.x*blockIdx.x +1;
	const int jj = blockIdx.y +1;
	const int kk = blockIdx.z +1;
	int nall = pitch.pDouble * (nvl[1]+2)*(nvl[2]+2);
	if((threadIdx.x == 0) && (blockIdx.x == 0) && (blockIdx.y == 0) && (blockIdx.z == 0))
	{
		/*printf("pitch: %d %d\n", (int)pitch.pDouble, (int)pitch.pInt);*/
		const int i = find_id(26,94,75);
		printf("device agent[tau] %g, %d, %d\n", (float)agent[tau*nall + i], i, nall);
	}
	if(ii < nvl[0]+1)
	{
        const int i = find_id(ii,jj,kk);
        if (!(type[i] & tissue(EMP)))
		{

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

#pragma unroll
		for(int s = -1; s <= 1; s+=2)
		{
			for (int d=0; d < 3; d+=1) {
			  const int j = find_id(ii +s*(d==0),jj +s*(d==1),kk +s*(d==2));

			  if (type[j] & tissue(EMP)) continue;

			  double del_phr = agent[phr * nall + i] - agent[phr * nall + j];

			  // diffusion of tau
			  double dum = 0.5 * (arr_prop.Dtau[ nall * d + i] + arr_prop.Dtau[nall * d + j]) * del_phr;
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
	}
	if((threadIdx.x == 0) && (blockIdx.x == 0) && (blockIdx.y == 0) && (blockIdx.z == 0))
	{
		/*printf("pitch: %d %d\n", (int)pitch.pDouble, (int)pitch.pInt);*/
		const int i = find_id(26,94,75);
		printf("device deriv[tau] %g, %d, %d\n", (float)deriv[tau*nall + i], i, nall);
	}
}

/* ----------------------------------------------------------------------*/
void ScenarioConnectomeStrategyCUDA::update() {

  using namespace ScenarioConnectomeAgents;

	static constexpr int BLOCK_DIM = 128;
	//const dim3 blocks(m_this->nvl[0]/BLOCK_DIM + (m_this->nvl[0]%BLOCK_DIM>0), m_this->nvl[1], m_this->nvl[2]);
	/*updateKernel<<<m_this->nall/BLOCK_DIM + (m_this->nall%BLOCK_DIM>0), BLOCK_DIM>>>(agent, deriv, type, arr_prop,m_this->dt, m_this->nall);*/
	const int nall = m_allocPitch.pDouble * (m_this->nvl[1]+2)*(m_this->nvl[2]+2);
	updateKernel<<<nall/BLOCK_DIM + (nall%BLOCK_DIM>0), BLOCK_DIM>>>(agent, deriv, type, arr_prop,m_this->dt, nall);

}

static __global__ void updateKernel(double* agent, const double* deriv, const int* type, const ScenarioConnectomeStrategyCUDA::array_properties arr_prop, double dt, int pnall)
{
	int i = threadIdx.x + blockDim.x*blockIdx.x;
  if(i < pnall) {
    if (type[i] & tissue(EMP)) return;
    
    // time integration (Euler's scheme)
    for (int ag_id=0; ag_id<ScenarioConnectomeAgents::num_agents; ag_id++)
      agent[ag_id * pnall + i] += deriv[ag_id * pnall + i] * dt;
  }
#if 0
	//const int jj = blockIdx.y +1;
	//const int kk = blockIdx.z +1;
	{
		int kk = i;
		int ii = kk%(nvl[0]+2);
		kk /= (nvl[0]+2);
		int jj = kk%(nvl[1]+2);
		kk /= (nvl[1]+2);
		if(ii == 26 && jj == 94 && kk==75)
		{
			const int a = find_id(ii,jj,kk);
			printf("device update agent[tau] %g %d (logic) %d, %d\n", (float)agent[tau*pitch.pDouble * (nvl[1]+2)*(nvl[2]+2) +a ], a, i, pnall);
		}
	}
  if(i < pnall) {
	i = (i/(nvl[0]+2))*pitch.pDouble + i%(nvl[0]+2);
	int nall = (pnall/(nvl[0]+2))*pitch.pDouble;

	{
		int kk = i;
		int ii = kk%pitch.pDouble;
		kk /= pitch.pDouble;
		int jj = kk%(nvl[1]+2);
		kk /= (nvl[1]+2);
		if(ii == 26 && jj == 94 && kk==75)
			printf("device update agent[tau] %g, %d, %d, %d\n", (float)agent[tau*nall + i], i, (int)(threadIdx.x + blockDim.x*blockIdx.x), nall);
	}
    
    if (type[i] & tissue(EMP)) return;
    
    // time integration (Euler's scheme)
    for (int ag_id=0; ag_id<ScenarioConnectomeAgents::num_agents; ag_id++)
      agent[ag_id * nall + i] += deriv[ag_id * nall + i] * dt;
  }
#endif
}
