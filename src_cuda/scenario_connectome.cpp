#include "scenario_connectome.h"

using namespace std;
using namespace ns_connectome;

/* ----------------------------------------------------------------------*/
ScenarioConnectome::ScenarioConnectome(int narg, char** arg, int rk, int np) {
  me = rk;
  nproc = np;

  MPI_Comm_split(MPI_COMM_WORLD,0,me,&world);

  allocations();

  if (!me) {
    printf("Reading input, setup the system ... \n");
    ofstream logfile;
    logfile.open (flog, ios::trunc);
    logfile << "*** LOG FILE ***" << endl;
    logfile.close();
  }

  scenario = arg[1];
  input->file(arg[2], this);

  // output initial step
  if (!me)
    printf("Writing output for the initial step ... \n");
  //output->lammpstrj(this);

  if (output->do_dump)
    output->dump(this);

  if (output->severy > 0)
    output->statistics(this);

  if (!me)
    printf("Integration started. \n");
  integrate(Nrun);

}

/* ----------------------------------------------------------------------*/
ScenarioConnectome::~ScenarioConnectome() {
  if(nim)
    nifti_image_free(nim);

  delete region;
  delete output;
  delete comm;
  delete init;
  delete input;
}

/* ----------------------------------------------------------------------*/
void ScenarioConnectome::allocations() {
  nvoxel = 0;
  nlocal = nghost = nall = 0;
  step = Nrun = 0;
  Nlog = 1000;

  dt = 0.0;
  nevery = -1;
  vlen = vlen_1 = vlen_2 = 0.0;
  vvol = vvol_1 = 0.0;

  for (int ag_id=0; ag_id<num_agents; ag_id++)
    init_val[ag_id] = -1.0;

  prop.Dtau_max = prop.diff_tau = 0.0;
  prop.dnt = 0.0;
  prop.D_sAb = prop.diff_sAb = 0.0;
  prop.D_mic = prop.diff_mic = 0.0;
  prop.cs = prop.sens_s = prop.cf = prop.sens_f = 0.0;
  prop.kp = prop.kn = 0.0;
  prop.ds = prop.df = 0.0;
  prop.es = 0.0;
  prop.Ha = 0.0;
  prop.ka = 0.0;

  prop.C_cir = 1.0;
  prop.c_cir = 0.0;
  prop.tau_cir = 1.0;
  prop.omega_cir = 0.0;

  prop.ktau = 0.0;
  prop.kphi = 0.0;
  prop.ephi = 0.0;

  input = new Input();
  init = new Init();
  comm = new Comm();
  output = new Output();
  region = new Region();

  nim = NULL;

  newton_flux = 1;

  // set tissue
  tissue.clear();
  tissue.resize(num_types);
  for (int i=0; i<num_types; i++)
    tissue[i] = 1 << i;
}

/* ----------------------------------------------------------------------*/
void ScenarioConnectome::integrate(int Nrun) {
  MPI_Barrier(world);
  double t0 = MPI_Wtime();
  double t1 = t0;

  double mu0,mu1, sig0,sig1;
  int N0,N1;

  mu0 = mu1 = sig0 = sig1 = 0.0;
  N0 = N1 = 0;

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
      if (!me) {
        double speed = float(Nlog)/(t2-t1);

        N1 = N0 + 1;
        mu1 = (N0 * mu0 + speed) / N1;
        sig1 = sqrt( N0 * sig0 * sig0 / N1
                   + N0 * (speed - mu0) * (speed - mu0) / (N1 * N1) );

        char buf[256];
        sprintf(buf,"Step %i, time lapsed = %g s, speed = %g +- %g steps/s \n",
                step, t2 - t0, mu1, sig1 );

        N0 = N1;
        mu0 = mu1;
        sig0 = sig1;

        printf(buf);

      }
      t1 = t2;
    }

  }

  MPI_Barrier(world);
  t1 = MPI_Wtime();
  if (!me) {
    char buf[256];
    sprintf(buf,"Final step: \n"
                "Total steps = %i \n"
                "Tot time (sec) = %g \n"
                "Speed (steps/s)= %g +- %g \n",
            Nrun, t1 - t0, mu1, sig1 );
    printf(buf);

    ofstream logfile;
    logfile.open (flog, ios::app);
    logfile << buf;

    logfile.close();
  }

}

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
          double dum = 0.5 * (prop.Dtau[d][i] + prop.Dtau[d][j]) * del_phr;
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
