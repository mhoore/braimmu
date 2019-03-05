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
  double t0,t1;

  MPI_Barrier(brn->world);
  t0 = MPI_Wtime();

  int iter = 0;
  while (iter < Nrun) {

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
      double t2;
      t2 = MPI_Wtime();
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
  int i,j,jj,ag_id;
  double del_sAb,del_fAb, del_mic;
  double dum;

  int nlocal = brn->nlocal;
  int nall = brn->nall;

  int *type = brn->type;

  double ***agent = brn->agent;

  int *num_neigh = brn->num_neigh;
  int **neigh = brn->neigh;

  // set derivatives of all voxels to zero
  for (ag_id=0; ag_id<num_agents; ag_id++)
    for (i=0; i<nall; i++)
      agent[ag_id][i][1] = 0.0;

  // spatial derivatives
  for (i=0; i<nlocal; i++) {
    if (type[i] == EMP_type) continue;

    for (jj=0; jj<num_neigh[i]; jj++) {
      j = neigh[i][jj];

      if (type[j] == EMP_type) continue;

      del_sAb = agent[sAb][i][0] - agent[sAb][j][0];

      // diffusion of sAb
      dum = brn->D_sAb * del_sAb;
      agent[sAb][i][1] -= dum;
      agent[sAb][j][1] += dum;

      // only in parenchyma
      if (type[i] != PAR_type) continue;
      if (type[j] != PAR_type) continue;

      del_fAb = agent[fAb][i][0] - agent[fAb][j][0];
      del_mic = agent[mic][i][0] - agent[mic][j][0];

      // migration of microglia toward higher sAb concentrations
      dum = brn->cs * del_sAb;
      if (del_sAb > 0.0)
        dum *= agent[mic][j][0];
      else
        dum *= agent[mic][i][0];

      agent[mic][i][1] += dum;
      agent[mic][j][1] -= dum;

      // migration of microglia toward higher fAb concentrations
      dum = brn->cf * del_fAb;
      if (del_fAb > 0.0)
        dum *= agent[mic][j][0];
      else
        dum *= agent[mic][i][0];

      agent[mic][i][1] += dum;
      agent[mic][j][1] -= dum;

      // diffusion of microglia
      dum = brn->D_mic * del_mic;
      agent[mic][i][1] -= dum;
      agent[mic][j][1] += dum;
    }

  }

}

/* ----------------------------------------------------------------------*/
void Run::update(Brain *brn) {
  int i,ag_id;
  double dum;

  int nlocal = brn->nlocal;

  int *type = brn->type;
  double dt = brn->dt;

  double ***agent = brn->agent;

  // update local voxels
  for (i=0; i<nlocal; i++) {
    if (type[i] == EMP_type) continue;

    dum = brn->kp * agent[sAb][i][0] * agent[fAb][i][0]
        + brn->kn * agent[sAb][i][0] * agent[sAb][i][0];

    // sAb
    agent[sAb][i][1] += agent[neu][i][0] * agent[cir][i][0]
                      - dum
                      - brn->ds * agent[mic][i][0] * agent[sAb][i][0];
    // fAb
    agent[fAb][i][1] += dum
                      - brn->df * agent[mic][i][0] * agent[fAb][i][0];

    // sAb & fAb efflux from CSF
    if (type[i] == CSF_type) {
      agent[sAb][i][1] -= brn->es * agent[sAb][i][0];
      agent[fAb][i][1] -= brn->es * agent[fAb][i][0];
    }

    // neuronal death due to astrogliosis
    agent[neu][i][1] -= (brn->dna * agent[ast][i][0]
                       + brn->dnf * agent[fAb][i][0]) * agent[neu][i][0];

    // astrogliosis
    dum = agent[fAb][i][0] * agent[mic][i][0];
    agent[ast][i][1] = brn->ka * (dum / (dum + brn->Ha) - agent[ast][i][0]);

    // time integration (Euler's scheme)
    for (ag_id=0; ag_id<num_agents; ag_id++)
      agent[ag_id][i][0] += agent[ag_id][i][1] * dt;

    //double rsq = sqrt(brn->x[i][0] * brn->x[i][0]
      //              + brn->x[i][1] * brn->x[i][1]
        //            + brn->x[i][2] * brn->x[i][2]);

    //if (rsq < 1.0e4)
      //agent[fAb][i][0] = 1.0e4 - rsq;

    // NOTE: correction: avoid negative values
    //if (agent[ag_id][i][0] < 0.0)
      //agent[ag_id][i][0] = 0.0;

    //printf("proc %i: HERE grad = %g, dir = %i tag= " TAGINT_FORMAT " \n",
    //       brn->me,grad[neu][j][dir],dir,brn->tag[i]);

  }

}
