#include "ScenarioConnectomeStrategyOACC.h"
#include "openacc.h"
#include "scenario_connectome.h"

using namespace std;

/* ----------------------------------------------------------------------*/
void ScenarioConnectomeStrategyOACC::derivatives() {

	using namespace ScenarioConnectomeAgents;

	auto &deriv = m_this->deriv;
	auto &type = m_this->type;
	auto &tissue = m_this->tissue;
	auto &nvl = m_this->nvl;
	auto &prop = m_this->prop;
	auto &agent = m_this->agent;
	const bool newton_flux = m_this->newton_flux;
	auto &arr_prop = m_this->arr_prop;

  // for (auto i = 0; i < num_agents; ++i) {
  //   m_this->agent_raw[i] = agent[i].data();
  //   m_this->deriv_raw[i] = deriv[i].data();
  // }

  // set derivatives of all voxels to zero
  for (int ag_id=0; ag_id<num_agents; ag_id++) {
    auto deriv_ptr = deriv[ag_id].data();
    auto agent_ptr = agent[ag_id].data();
    const auto N = deriv[ag_id].size();

    #pragma acc enter data create(deriv_ptr[0:N], agent_ptr[0:N])
    #pragma acc parallel loop present(deriv_ptr, agent_ptr)
    for (auto i = 0; i < N; ++i){
      deriv_ptr[i] = 0.0;
      agent_ptr[i] = 0.0;
    }
  }

  auto deriv_sAb = deriv[sAb].data();
  auto deriv_fAb = deriv[fAb].data();
  auto deriv_phr = deriv[phr].data();
  auto deriv_tau = deriv[tau].data();
  auto deriv_neu = deriv[neu].data();
  auto deriv_ast = deriv[ast].data();
  auto deriv_cir = deriv[cir].data();
  auto deriv_mic = deriv[mic].data();

  auto agent_sAb = agent[sAb].data();
  auto agent_fAb = agent[fAb].data();
  auto agent_phr = agent[phr].data();
  auto agent_tau = agent[tau].data();
  auto agent_neu = agent[neu].data();
  auto agent_ast = agent[ast].data();
  auto agent_cir = agent[cir].data();
  auto agent_mic = agent[mic].data();

  // spatial derivatives
  for (int p=0; p<2; p++)
    { //printf("%d",omp_get_num_threads());
    //#pragma omp for collapse(2)

      int k_dim = nvl[2] + 1;
      int j_dim = nvl[1] + 1;
    #pragma acc parallel loop collapse(2) present(deriv_sAb, deriv_fAb, deriv_phr, deriv_tau, deriv_neu, deriv_ast, deriv_cir, deriv_mic, agent_sAb, agent_fAb, agent_phr, agent_tau, agent_neu, agent_ast, agent_cir, agent_mic)
    for (int kk=1; kk<k_dim; kk++)
      for (int jj=1; jj<j_dim; jj++)
        //f r (int ii=1; ii<nvl[0]+1; ii+2) {
        for (int ii = (kk^jj^p)&1; ii<nvl[0]+1; ii+=2){
          //cin >> ii;
          int i = m_this->find_id(ii,jj,kk);
          if (type[i] & tissue[EMP]) continue;

          // direct function or time derivatives

          // sAb, fAb, and tau efflux from CSF
          if (type[i] & tissue[CSF]) {
            deriv_sAb[i] -= prop.es * agent_sAb[i];
            deriv_fAb[i] -= prop.es * agent_fAb[i];
            deriv_phr[i] -= prop.ephi * agent_phr[i];
          }
          // in parenchyma (WM and GM)
          else {
            double dum = prop.kp * agent_sAb[i] * agent_fAb[i]
                       + prop.kn * agent_sAb[i] * agent_sAb[i];

            // sAb
            deriv_sAb[i] += agent_neu[i] * agent_cir[i]
                              - dum
                              - prop.ds * agent_mic[i] * agent_sAb[i];
            // fAb
            deriv_fAb[i] += dum
                           - prop.df * agent_mic[i] * agent_fAb[i];

            dum = prop.ktau * agent_phr[i];

            // tau protein phosphorylation due to fAb and neu
            deriv_phr[i] += prop.kphi * agent_fAb[i] * agent_neu[i]
                           - dum;

            // tau tangle formation from phosphorylated tau
            deriv[tau][i] += dum;

            // neuronal death due to tau aggregation
            deriv_neu[i] -= prop.dnt * agent_tau[i] * agent_neu[i];

            // astrogliosis
            dum = agent_fAb[i] * agent_mic[i];
            deriv_ast[i] = prop.ka * (dum / (dum + prop.Ha) - agent_ast[i]);

            // circadian rhythm
            if (prop.c_cir > 0)
              deriv_cir[i] = - prop.C_cir * prop.c_cir * prop.omega_cir
                              * sin(prop.omega_cir * m_this->dt * m_this->step);
          }

          // spatial derivatives: fluxes
          int n_ngh, ngh[6];
          ngh[0] = m_this->find_id(ii+1,jj,kk);
          ngh[1] = m_this->find_id(ii,jj+1,kk);
          ngh[2] = m_this->find_id(ii,jj,kk+1);
          n_ngh = 3;

          /*
          if (!newton_flux) {
            ngh[3] = m_this->find_id(ii-1,jj,kk);
            ngh[4] = m_this->find_id(ii,jj-1,kk);
            ngh[5] = m_this->find_id(ii,jj,kk-1);
            n_ngh = 6;
          }
          */
          for (int c=0; c<n_ngh; ++c) {
            int j = ngh[c];
            int d = c;
           // if (c >= 3)
              //d = c - 3;

            if (type[j] & tissue[EMP]) continue;

            double del_phr = agent_phr[i] - agent_phr[j];

            // diffusion of tau
            double dum = 0.5 * (arr_prop.Dtau[d][i] + arr_prop.Dtau[d][j]) * del_phr;
            deriv_phr[i] -= dum;
            // if (newton_flux)
            deriv_phr[j] += dum;

            double del_sAb = agent_sAb[i] - agent_sAb[j];

            // diffusion of sAb
            dum = prop.D_sAb * del_sAb;
            deriv_sAb[i] -= dum;
            // if (newton_flux)
            deriv_sAb[j] += dum;

            // only in parenchyma
            if (type[i] & tissue[WM] || type[i] & tissue[GM])
              if (type[j] & tissue[WM] || type[j] & tissue[GM]) {
                double del_fAb = agent_fAb[i] - agent_fAb[j];
                double del_mic = agent_mic[i] - agent_mic[j];

                // migration of microglia toward higher sAb concentrations
                dum = prop.cs * del_sAb;
                if (del_sAb > 0.0)
                  dum *= agent_mic[j];
                else
                  dum *= agent_mic[i];

                deriv_mic[i] += dum;
                //if (newton_flux)
                deriv_mic[j] -= dum;

                // migration of microglia toward higher fAb concentrations
                dum = prop.cf * del_fAb;
                if (del_fAb > 0.0)
                  dum *= agent_mic[j];
                else
                  dum *= agent_mic[i];

                deriv_mic[i] += dum;
                //if (newton_flux)
                deriv_mic[j] -= dum;

                // diffusion of microglia
                dum = prop.D_mic * del_mic;
                deriv_mic[i] -= dum;
                //if (newton_flux)
                deriv_mic[j] += dum;
              }
          }
    /*
    #pragma omp for collapse(3)
    for (int kk=1; kk<nvl[2]+1; kk++)
      for (int jj=1; jj<nvl[1]+1; jj++)
         for (int ii=1; ii<nvl[0]+1; ii++){
             deriv_phr[i]= 2buff_deriv_phr[0][i]+ 2buff_deriv_phr[1][i] 
             deriv_sAb[i]= 2buff_deriv_sAb[0][i]+ 2buff_deriv_sAb[1][i]
             deriv_mic[i]= 2buff_deriv_mic[0][i]+ 2buff_deriv_mic[1][i]
             }*/
        }
    }

  for (int ag_id=0; ag_id<num_agents; ag_id++) {
    auto deriv_ptr = deriv[ag_id].data();
    auto agent_ptr = agent[ag_id].data();
    const auto N = deriv[ag_id].size();
    #pragma acc exit data copyout(deriv_ptr[0:N], agent_ptr[0:N])
  }

}

/* ----------------------------------------------------------------------*/
void ScenarioConnectomeStrategyOACC::update() {

	using namespace ScenarioConnectomeAgents;
  //#pragma omp parallel
  //{
  // update local voxels
  // #pragma acc parallel loop collapse(2)
  for (int ag_id=0; ag_id<num_agents; ag_id++) 
    for (int kk=1; kk<m_this->nvl[2]+1; kk++)
      for (int jj=1; jj<m_this->nvl[1]+1; jj++)
        for (int ii=1; ii<m_this->nvl[0]+1; ii++) {
          int i = m_this->find_id(ii,jj,kk);
          if (m_this->type[i] & m_this->tissue[EMP]) continue;

          // time integration (Euler's scheme)
           m_this->agent[ag_id][i] += m_this->deriv[ag_id][i] * m_this->dt;
      }

  //}
}
