#include "scenario_connectome.h"

using namespace std;

/* ----------------------------------------------------------------------*/
void ScenarioConnectome::derivatives() {

  // set derivatives of all voxels to zero
  for (int ag_id=0; ag_id<num_agents; ag_id++)
    fill(deriv[ag_id].begin(), deriv[ag_id].end(), 0.);

  // spatial derivatives
  for (int kk=1; kk<nvl[2]+1; kk++)
    for (int jj=1; jj<nvl[1]+1; jj++)
      for (int ii=1; ii<nvl[0]+1; ii++) {
        int i = find_id(ii,jj,kk);
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
                            * sin(prop.omega_cir * dt * step);
        }

        // spatial derivatives: fluxes
        int n_ngh, ngh[6];
        ngh[0] = find_id(ii+1,jj,kk);
        ngh[1] = find_id(ii,jj+1,kk);
        ngh[2] = find_id(ii,jj,kk+1);
        n_ngh = 3;

        if (!newton_flux) {
          ngh[3] = find_id(ii-1,jj,kk);
          ngh[4] = find_id(ii,jj-1,kk);
          ngh[5] = find_id(ii,jj,kk-1);
          n_ngh = 6;
        }

        for (int c=0; c<n_ngh; ++c) {
          int j = ngh[c];
          int d = c;
          if (c >= 3)
            d = c - 3;

          if (type[j] & tissue[EMP]) continue;

          double del_phr = agent[phr][i] - agent[phr][j];

          // diffusion of tau
          double dum = 0.5 * (arr_prop.Dtau[d][i] + arr_prop.Dtau[d][j]) * del_phr;
          deriv[phr][i] -= dum;
          if (newton_flux)
            deriv[phr][j] += dum;

          double del_sAb = agent[sAb][i] - agent[sAb][j];

          // diffusion of sAb
          dum = prop.D_sAb * del_sAb;
          deriv[sAb][i] -= dum;
          if (newton_flux)
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
              if (newton_flux)
                deriv[mic][j] -= dum;

              // migration of microglia toward higher fAb concentrations
              dum = prop.cf * del_fAb;
              if (del_fAb > 0.0)
                dum *= agent[mic][j];
              else
                dum *= agent[mic][i];

              deriv[mic][i] += dum;
              if (newton_flux)
                deriv[mic][j] -= dum;

              // diffusion of microglia
              dum = prop.D_mic * del_mic;
              deriv[mic][i] -= dum;
              if (newton_flux)
                deriv[mic][j] += dum;
            }
        }
      }

}

/* ----------------------------------------------------------------------*/
void ScenarioConnectome::update() {

  // update local voxels
  for (int kk=1; kk<nvl[2]+1; kk++)
    for (int jj=1; jj<nvl[1]+1; jj++)
      for (int ii=1; ii<nvl[0]+1; ii++) {
        int i = find_id(ii,jj,kk);
        if (type[i] & tissue[EMP]) continue;

        // time integration (Euler's scheme)
        for (int ag_id=0; ag_id<num_agents; ag_id++)
          agent[ag_id][i] += deriv[ag_id][i] * dt;
      }

}
