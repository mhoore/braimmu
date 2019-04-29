#include <mpi.h>
#include "math.h"

#include "pointers.h"
#include "brain.h"
#include "run.h"

using namespace std;
using namespace brain_NS;

/* ---------------------------------------------------------------------- */
Run::Run() {
}

/* ----------------------------------------------------------------------*/
Run::~Run() {
}

/* ----------------------------------------------------------------------*/
void Run::integrate(Brain *brn, int Nrun) {
  MPI_Barrier(brn->world);
  double t0 = MPI_Wtime();
  double t1 = t0;

  int iter = 0;
  while (iter < Nrun) {

    brn->comm->forward_comm(brn);

    derivatives(brn);

    brn->comm->reverse_comm(brn);

    // communicate bond connections

    update(brn);
    brn->step++;
    iter++;

    if (brn->output->do_dump)
      brn->output->dump(brn);
    if (brn->output->severy > 0)
      brn->output->statistics(brn);

    if (brn->step % brn->Nlog == 0) {
      MPI_Barrier(brn->world);
      double t2 = MPI_Wtime();
      if (!brn->me)
        printf("Step %i, time lapsed = %g sec, speed = %g steps/sec \n",
               brn->step, t2 - t0, float(brn->Nlog)/(t2-t1) );
      t1 = t2;
    }

  }

  MPI_Barrier(brn->world);
  t1 = MPI_Wtime();
  if (!brn->me)
    printf("Final step, total iterations = %i, total time lapsed = %g sec, speed = %g steps/sec \n",
           Nrun, t1 - t0, float(Nrun)/(t1-t0) );

}

/* ----------------------------------------------------------------------*/
void Run::derivatives(Brain *brn) {
  int nall = brn->nall;

  int *type = brn->type;

  double **agent = brn->agent;
  double **deriv = brn->deriv;

  int *num_neigh = brn->num_neigh;
  int **neigh = brn->neigh;

  int *nvl = brn->nvl;

  // set derivatives of all voxels to zero
  for (int ag_id=0; ag_id<num_agents; ag_id++)
    for (int i=0; i<nall; i++)
      deriv[ag_id][i] = 0.0;

  // spatial derivatives
  for (int kk=1; kk<nvl[2]+1; kk++)
    for (int jj=1; jj<nvl[1]+1; jj++)
      for (int ii=1; ii<nvl[0]+1; ii++) {
        int i = find_id(brn,ii,jj,kk);
        if (type[i] == EMP_type) continue;

        for (int hh=0; hh<num_neigh[i]; hh++) {
          int j = neigh[i][hh];

          if (type[j] == EMP_type) continue;

          double del_sAb = agent[sAb][i] - agent[sAb][j];

          // diffusion of sAb
          double dum = brn->D_sAb * del_sAb;
          deriv[sAb][i] -= dum;
          deriv[sAb][j] += dum;

          // only in parenchyma
          if (type[i] != WM_type && type[i] != GM_type) continue;
          if (type[j] != WM_type && type[j] != GM_type) continue;

          double del_fAb = agent[fAb][i] - agent[fAb][j];
          double del_mic = agent[mic][i] - agent[mic][j];

          // migration of microglia toward higher sAb concentrations
          dum = brn->cs * del_sAb;
          if (del_sAb > 0.0)
            dum *= agent[mic][j];
          else
            dum *= agent[mic][i];

          deriv[mic][i] += dum;
          deriv[mic][j] -= dum;

          // migration of microglia toward higher fAb concentrations
          dum = brn->cf * del_fAb;
          if (del_fAb > 0.0)
            dum *= agent[mic][j];
          else
            dum *= agent[mic][i];

          deriv[mic][i] += dum;
          deriv[mic][j] -= dum;

          // diffusion of microglia
          dum = brn->D_mic * del_mic;
          deriv[mic][i] -= dum;
          deriv[mic][j] += dum;
        }
      }

}

/* ----------------------------------------------------------------------*/
void Run::update(Brain *brn) {
  int *type = brn->type;
  double dt = brn->dt;

  double **agent = brn->agent;
  double **deriv = brn->deriv;

  int *nvl = brn->nvl;

  // update local voxels
  for (int kk=1; kk<nvl[2]+1; kk++)
    for (int jj=1; jj<nvl[1]+1; jj++)
      for (int ii=1; ii<nvl[0]+1; ii++) {
        int i = find_id(brn,ii,jj,kk);
        if (type[i] == EMP_type) continue;

        double dum = brn->kp * agent[sAb][i] * agent[fAb][i]
                   + brn->kn * agent[sAb][i] * agent[sAb][i];

        // sAb
        deriv[sAb][i] += agent[neu][i] * agent[cir][i]
                          - dum
                          - brn->ds * agent[mic][i] * agent[sAb][i];
        // fAb
        deriv[fAb][i] += dum
                          - brn->df * agent[mic][i] * agent[fAb][i];

        // sAb & fAb efflux from CSF
        if (type[i] == CSF_type) {
          deriv[sAb][i] -= brn->es * agent[sAb][i];
          deriv[fAb][i] -= brn->es * agent[fAb][i];
        }

        // neuronal death due to astrogliosis
        deriv[neu][i] -= (brn->dna * agent[ast][i]
                           + brn->dnf * agent[fAb][i]) * agent[neu][i];

        // astrogliosis
        dum = agent[fAb][i] * agent[mic][i];
        deriv[ast][i] = brn->ka * (dum / (dum + brn->Ha) - agent[ast][i]);

        // circadian rhythm
        if (brn->c_cir > 0)
          deriv[cir][i] = - brn->C_cir * brn->c_cir * brn->omega_cir
                            * sin(brn->omega_cir * dt * brn->step);

        // time integration (Euler's scheme)
        for (int ag_id=0; ag_id<num_agents; ag_id++)
          agent[ag_id][i] += deriv[ag_id][i] * dt;

      }

}

/* ----------------------------------------------------------------------
 * Find the local voxel id from local coordinates i,j,k
 * ----------------------------------------------------------------------*/
int Run::find_id(Brain *brn, int i, int j, int k) {
  int *nvl = brn->nvl;
  return i + (nvl[0] + 2) * (j + (nvl[1] + 2) * k);
}

/*    ////////DEBUG/////////////////////
    FILE* fw;
    fw = fopen("out.txt","a");
    for (int i=0; i<brn->nall; i++) {
      fprintf(fw,"proc%i: %i %i %g %g %g ",
              brn->me, brn->step, i, brn->x[i][0], brn->x[i][1], brn->x[i][2]);
      for (int ag_id=0; ag_id<num_agents; ag_id++)
        fprintf(fw,"%s %g %g ",
                ag_str[ag_id].c_str(), brn->agent[ag_id][i][0], brn->agent[ag_id][i][1]);
      fprintf(fw,"\n");
    }
    fclose(fw);
    ////////DEBUG/////////////////////
*/
