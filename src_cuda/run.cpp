#include <mpi.h>
#include "math.h"

#include "pointers.h"
#include "brain.h"

using namespace std;
using namespace brain_NS;

/* ----------------------------------------------------------------------*/
void Brain::integrate(int Nrun) {
  MPI_Barrier(world);
  double t0 = MPI_Wtime();
  double t1 = t0;

  int iter = 0;
  while (iter < Nrun) {

    comm->forward_comm(this);

    derivatives();

    if (newton_flux)
      comm->reverse_comm(this);

    // communicate bond connections

    update();
    step++;
    iter++;

    if (output->do_dump)
      output->dump(this);
    if (output->severy > 0)
      output->statistics(this);

    if (step % Nlog == 0) {
      MPI_Barrier(world);
      double t2 = MPI_Wtime();
      if (!me)
        printf("Step %i, time lapsed = %g sec, speed = %g steps/sec \n",
               step, t2 - t0, float(Nlog)/(t2-t1) );
      t1 = t2;
    }

  }

  MPI_Barrier(world);
  t1 = MPI_Wtime();
  if (!me)
    printf("Final step, total iterations = %i, total time lapsed = %g sec, speed = %g steps/sec \n",
           Nrun, t1 - t0, float(Nrun)/(t1-t0) );

}

/* ----------------------------------------------------------------------*/
void Brain::derivatives() {

  // set derivatives of all voxels to zero
  for (int ag_id=0; ag_id<num_agents; ag_id++)
    for (int i=0; i<nall; i++)
      deriv[ag_id][i] = 0.0;

  // spatial derivatives
  for (int kk=1; kk<nvl[2]+1; kk++)
    for (int jj=1; jj<nvl[1]+1; jj++)
      for (int ii=1; ii<nvl[0]+1; ii++) {
        int i = find_id(ii,jj,kk);
        if (type[i] == EMP_type) continue;

        double dum = kp * agent[sAb][i] * agent[fAb][i]
                   + kn * agent[sAb][i] * agent[sAb][i];

        // sAb
        deriv[sAb][i] += agent[neu][i] * agent[cir][i]
                          - dum
                          - ds * agent[mic][i] * agent[sAb][i];
        // fAb
        deriv[fAb][i] += dum
                          - df * agent[mic][i] * agent[fAb][i];

        // sAb & fAb efflux from CSF
        if (type[i] == CSF_type) {
          deriv[sAb][i] -= es * agent[sAb][i];
          deriv[fAb][i] -= es * agent[fAb][i];
        }

        // neuronal death due to astrogliosis
        deriv[neu][i] -= (dna * agent[ast][i]
                           + dnf * agent[fAb][i]) * agent[neu][i];

        // astrogliosis
        dum = agent[fAb][i] * agent[mic][i];
        deriv[ast][i] = ka * (dum / (dum + Ha) - agent[ast][i]);

        /// fluxes
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

          if (type[j] == EMP_type) continue;

          double del_sAb = agent[sAb][i] - agent[sAb][j];

          // diffusion of sAb
          double dum = D_sAb * del_sAb;
          deriv[sAb][i] -= dum;
          if (newton_flux)
            deriv[sAb][j] += dum;

          // only in parenchyma
          if (type[i] != WM_type && type[i] != GM_type) continue;
          if (type[j] != WM_type && type[j] != GM_type) continue;

          double del_fAb = agent[fAb][i] - agent[fAb][j];
          double del_mic = agent[mic][i] - agent[mic][j];

          // migration of microglia toward higher sAb concentrations
          dum = cs * del_sAb;
          if (del_sAb > 0.0)
            dum *= agent[mic][j];
          else
            dum *= agent[mic][i];

          deriv[mic][i] += dum;
          if (newton_flux)
            deriv[mic][j] -= dum;

          // migration of microglia toward higher fAb concentrations
          dum = cf * del_fAb;
          if (del_fAb > 0.0)
            dum *= agent[mic][j];
          else
            dum *= agent[mic][i];

          deriv[mic][i] += dum;
          if (newton_flux)
            deriv[mic][j] -= dum;

          // diffusion of microglia
          dum = D_mic * del_mic;
          deriv[mic][i] -= dum;
          if (newton_flux)
            deriv[mic][j] += dum;
        }
      }

}

/* ----------------------------------------------------------------------*/
void Brain::update() {

  // update local voxels
  for (int kk=1; kk<nvl[2]+1; kk++)
    for (int jj=1; jj<nvl[1]+1; jj++)
      for (int ii=1; ii<nvl[0]+1; ii++) {
        int i = find_id(ii,jj,kk);
        if (type[i] == EMP_type) continue;

        // time integration (Euler's scheme)
        for (int ag_id=0; ag_id<num_agents; ag_id++)
          agent[ag_id][i] += deriv[ag_id][i] * dt;

      }

}

/* ----------------------------------------------------------------------
 * Find the local voxel id from local coordinates i,j,k
 * ----------------------------------------------------------------------*/
int Brain::find_id(int i, int j, int k) {
  return i + (nvl[0] + 2) * (j + (nvl[1] + 2) * k);
}
