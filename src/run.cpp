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

        // direct function or time derivatives

        // sAb, fAb, and tau efflux from CSF
        if (type[i] == CSF_type) {
          deriv[sAb][i] -= es * agent[sAb][i];
          deriv[fAb][i] -= es * agent[fAb][i];
          deriv[phr][i] -= ephi * agent[phr][i];
        }
        // in parenchyma (WM and GM)
        else {
          double dum = kp * agent[sAb][i] * agent[fAb][i]
                     + kn * agent[sAb][i] * agent[sAb][i];

          // sAb
          deriv[sAb][i] += agent[neu][i] * agent[cir][i]
                            - dum
                            - ds * agent[mic][i] * agent[sAb][i];
          // fAb
          deriv[fAb][i] += dum
                           - df * agent[mic][i] * agent[fAb][i];

          dum = ktau * agent[phr][i];

          // tau protein phosphorylation due to fAb and neu
          deriv[phr][i] += kphi * agent[fAb][i] * agent[neu][i]
                         - dum;

          // tau tangle formation from phosphorylated tau
          deriv[tau][i] += dum;

          // neuronal death due to tau aggregation
          deriv[neu][i] -= dnt * agent[tau][i] * agent[neu][i];

          // astrogliosis
          dum = agent[fAb][i] * agent[mic][i];
          deriv[ast][i] = ka * (dum / (dum + Ha) - agent[ast][i]);

          // circadian rhythm
          if (c_cir > 0)
            deriv[cir][i] = - C_cir * c_cir * omega_cir
                            * sin(omega_cir * dt * step);
        }

        // spatial derivatives: fluxes
        for (int hh=0; hh<num_neigh[i]; hh++) {
          int j = neigh[i][hh];

          if (type[j] == EMP_type) continue;

          double del_phr = agent[phr][i] - agent[phr][j];

          // diffusion of tau
          double dum = 0.5 * (Dtau[hh][i] + Dtau[hh][j]) * del_phr;
          deriv[phr][i] -= dum;
          deriv[phr][j] += dum;

          double del_sAb = agent[sAb][i] - agent[sAb][j];

          // diffusion of sAb
          dum = D_sAb * del_sAb;
          deriv[sAb][i] -= dum;
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
          deriv[mic][j] -= dum;

          // migration of microglia toward higher fAb concentrations
          dum = cf * del_fAb;
          if (del_fAb > 0.0)
            dum *= agent[mic][j];
          else
            dum *= agent[mic][i];

          deriv[mic][i] += dum;
          deriv[mic][j] -= dum;

          // diffusion of microglia
          dum = D_mic * del_mic;
          deriv[mic][i] -= dum;
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

/*    ////////DEBUG/////////////////////
    FILE* fw;
    fw = fopen("out.txt","a");
    for (int i=0; i<nall; i++) {
      fprintf(fw,"proc%i: %i %i %g %g %g ",
              me, step, i, x[i][0], x[i][1], x[i][2]);
      for (int ag_id=0; ag_id<num_agents; ag_id++)
        fprintf(fw,"%s %g %g ",
                ag_str[ag_id].c_str(), agent[ag_id][i][0], agent[ag_id][i][1]);
      fprintf(fw,"\n");
    }
    fclose(fw);
    ////////DEBUG/////////////////////
*/
