#include <mpi.h>
#include "math.h"

#include "pointers.h"
#include "brain.h"

using namespace std;
using namespace brain_NS;

/* ----------------------------------------------------------------------*/
void Brain::integrate(int Nrun) {
  double t0,t1;

  MPI_Barrier(world);
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
      double t2;
      t2 = MPI_Wtime();
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
  int i,j,jj,ag_id;
  double del_sAb,del_fAb, del_mic;
  double dum;

  // set derivatives of all voxels to zero
  for (ag_id=0; ag_id<num_agents; ag_id++)
    for (i=0; i<nall; i++)
      deriv[ag_id][i] = 0.0;

  // spatial derivatives
  for (i=0; i<nlocal; i++) {
    if (type[i] == EMP_type) continue;

    //int ii,jj,kk;
    //tagint itag = brn->init->find_me(brn,ii,jj,kk);
    //j = brn->map[itag];

    for (jj=0; jj<num_neigh[i]; jj++) {
      j = neigh[i][jj];

      if (type[j] == EMP_type) continue;

      del_sAb = agent[sAb][i] - agent[sAb][j];

      // diffusion of sAb
      dum = D_sAb * del_sAb;
      deriv[sAb][i] -= dum;
      deriv[sAb][j] += dum;

      // only in parenchyma
      if (type[i] != WM_type && type[i] != GM_type) continue;
      if (type[j] != WM_type && type[j] != GM_type) continue;

      del_fAb = agent[fAb][i] - agent[fAb][j];
      del_mic = agent[mic][i] - agent[mic][j];

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
  int i,ag_id;
  double dum;

  // update local voxels
  for (i=0; i<nlocal; i++) {
    if (type[i] == EMP_type) continue;

    dum = kp * agent[sAb][i] * agent[fAb][i]
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

    // time integration (Euler's scheme)
    for (ag_id=0; ag_id<num_agents; ag_id++)
      agent[ag_id][i] += deriv[ag_id][i] * dt;

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
