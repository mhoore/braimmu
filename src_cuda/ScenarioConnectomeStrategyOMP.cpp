#include "ScenarioConnectomeStrategyOMP.h"
#include "omp.h"
#include "scenario_connectome.h"

using namespace std;

/* ----------------------------------------------------------------------*/
void ScenarioConnectomeStrategyOMP::derivatives() {

	using namespace ScenarioConnectomeAgents;

	auto &deriv = m_this->deriv;
        //auto &2buff_deriv = m_this->2buff_deriv;
	auto &type = m_this->type;
	auto &tissue = m_this->tissue;
	auto &nvl = m_this->nvl;
	auto &prop = m_this->prop;
	auto &agent = m_this->agent;
	const bool newton_flux = m_this->newton_flux;
	auto &arr_prop = m_this->arr_prop;

  // set derivatives of all voxels to zero
  //#pragma omp parallel
  #pragma omp parallel
  {
  for (int ag_id=0; ag_id<num_agents; ag_id++)
    fill(deriv[ag_id].begin(), deriv[ag_id].end(), 0.);
   // fill(2buff_deriv[ag_id].begin(), odd_deriv[ag_id].end(), 0.);
  }
  // spatial derivatives
  for (int p=0; p<2; p++)
    #pragma omp parallel
    { //printf("%d",omp_get_num_threads());
    #pragma omp for collapse(2)
    for (int kk=1; kk<nvl[2]+1; kk++)
      for (int jj=1; jj<nvl[1]+1; jj++)
        //f r (int ii=1; ii<nvl[0]+1; ii+2) {
        for (int ii = (kk^jj^p)&1; ii<nvl[0]+1; ii+=2){
          //cin >> ii;
          int i = m_this->find_id(ii,jj,kk);
          if (type[i] & tissue[EMP]) continue;

          // direct function or time derivatives

          // sAb, fAb, and tau efflux from CSF
          if (type[i] & tissue[CSF]) {
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

            double del_phr = agent[phr][i] - agent[phr][j];

            // diffusion of tau
            double dum = 0.5 * (arr_prop.Dtau[d][i] + arr_prop.Dtau[d][j]) * del_phr;
            deriv[phr][i] -= dum;
            // if (newton_flux)
            deriv[phr][j] += dum;

            double del_sAb = agent[sAb][i] - agent[sAb][j];

            // diffusion of sAb
            dum = prop.D_sAb * del_sAb;
            deriv[sAb][i] -= dum;
            // if (newton_flux)
            deriv[sAb][j] += dum;

            // only in parenchyma
            if (type[i] & tissue[WM] || type[i] & tissue[GM])
              if (type[j] & tissue[WM] || type[j] & tissue[GM]) {
                double del_fAb = agent[fAb][i] - agent[fAb][j];
                double del_mic = agent[mic][i] - agent[mic][j];

                // migration of microglia toward higher sAb concentrations
                dum = prop.cs * del_sAb;
                if (del_sAb > 0.0)
                  dum *= agent[mic][j];
                else
                  dum *= agent[mic][i];

                deriv[mic][i] += dum;
                //if (newton_flux)
                deriv[mic][j] -= dum;

                // migration of microglia toward higher fAb concentrations
                dum = prop.cf * del_fAb;
                if (del_fAb > 0.0)
                  dum *= agent[mic][j];
                else
                  dum *= agent[mic][i];

                deriv[mic][i] += dum;
                //if (newton_flux)
                deriv[mic][j] -= dum;

                // diffusion of microglia
                dum = prop.D_mic * del_mic;
                deriv[mic][i] -= dum;
                //if (newton_flux)
                deriv[mic][j] += dum;
              }
          }
    /*
    #pragma omp for collapse(3)
    for (int kk=1; kk<nvl[2]+1; kk++)
      for (int jj=1; jj<nvl[1]+1; jj++)
         for (int ii=1; ii<nvl[0]+1; ii++){
             deriv[phr][i]= 2buff_deriv[phr][0][i]+ 2buff_deriv[phr][1][i] 
             deriv[sAb][i]= 2buff_deriv[sAb][0][i]+ 2buff_deriv[sAb][1][i]
             deriv[mic][i]= 2buff_deriv[mic][0][i]+ 2buff_deriv[mic][1][i]
             }*/
           }
        }
  }

/* ----------------------------------------------------------------------*/
void ScenarioConnectomeStrategyOMP::update() {

	using namespace ScenarioConnectomeAgents;
  #pragma omp parallel
  {
  // update local voxels
  #pragma omp for collapse(2)
  for (int ag_id=0; ag_id<num_agents; ag_id++) 
    for (int kk=1; kk<m_this->nvl[2]+1; kk++)
      for (int jj=1; jj<m_this->nvl[1]+1; jj++)
        for (int ii=1; ii<m_this->nvl[0]+1; ii++) {
          int i = m_this->find_id(ii,jj,kk);
          if (m_this->type[i] & m_this->tissue[EMP]) continue;

          // time integration (Euler's scheme)
           m_this->agent[ag_id][i] += m_this->deriv[ag_id][i] * m_this->dt;
      }

  }
}
