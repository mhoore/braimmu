#include <mpi.h>
#include "math.h"

#include "pointers.h"
#include "brain.h"
#include "input.h"
#include "output.h"

using namespace std;
using namespace brain_NS;

/* ----------------------------------------------------------------------*/
Brain::Brain(int narg, char **arg, int rk, int np) {
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
  input->file(arg[1], this);

  // output initial step
  if (!me)
    printf("Writing output for the initial step ... \n");
  output->lammpstrj(this);

  if (output->do_dump)
    output->dump(this);

  if (output->severy > 0)
    output->statistics(this);

  if (!me)
    printf("Integration started. \n");
  integrate(Nrun);

}

/* ----------------------------------------------------------------------*/
Brain::~Brain() {
  if(nim)
    nifti_image_free(nim);

  delete region;
  delete output;
  delete comm;
  delete init;
  delete input;
}

/* ----------------------------------------------------------------------*/
void Brain::allocations() {
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

  Dtau_max = diff_tau = 0.0;
  ktau = ephi = kphi = dnt = 0.0;
  D_sAb = diff_sAb = 0.0;
  D_mic = diff_mic = 0.0;
  cs = sens_s = cf = sens_f = 0.0;
  kp = kn = 0.0;
  ds = df = 0.0;
  es = 0.0;

  C_cir = 1.0;
  c_cir = 0.0;
  tau_cir = 1.0;

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
