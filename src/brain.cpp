#include <mpi.h>
#include "math.h"

#include "pointers.h"
#include "brain.h"
#include "input.h"
#include "memory.h"
#include "output.h"

using namespace std;
using namespace brain_NS;

/* ----------------------------------------------------------------------*/
Brain::Brain(int narg, char **arg, int rk, int np) {
  me = rk;
  nproc = np;

  MPI_Comm_split(MPI_COMM_WORLD,0,me,&world);

  allocations();

  if (!me)
    printf("Reading input, setup the system ... \n");
  input->file(arg[1], this);

  // output initial step
  if (!me)
    printf("Writing output for the initial step ... \n");
//  output->lammpstrj(this);

  if (output->do_dump)
    output->dump(this);

  if (output->severy > 0)
    output->statistics(this);

  if (!me)
    printf("Integration started. \n");
  run->integrate(this, Nrun);

  //printf("proc %i: xlo = %g \n", me, xlo);
  //MPI_Barrier(MPI_COMM_WORLD);
  //printf("proc %i here3 \n",brn->me);
  //if (brn->me == 2)
    //printf("proc %i: unpack itag = %li, c=%i \n",brn->me,itag,c-1);

}

/* ----------------------------------------------------------------------*/
Brain::~Brain() {
  destroy();

  delete region;
  delete output;
  delete comm;
  delete run;
  delete init;
  delete input;
  delete memory;
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

  D_sAb = diff_sAb = 0.0;
  D_mic = diff_mic = 0.0;
  cs = sens_s = cf = sens_f = 0.0;
  kp = kn = 0.0;
  ds = df = 0.0;
  es = 0.0;

  C_cir = 1.0;
  c_cir = 0.0;
  tau_cir = 1.0;

  memory->create(nv,3,"nv");
  memory->create(nvl,3,"nvl");
  memory->create(npart,3,"npart");

  memory->create(boxlo,3,"boxlo");
  memory->create(boxhi,3,"boxhi");
  memory->create(lbox,3,"lbox");

  memory->create(xlo,3,"xlo");
  memory->create(xhi,3,"xhi");

  memory = new Memory();
  input = new Input();
  init = new Init();
  run = new Run();
  comm = new Comm();
  output = new Output();
  region = new Region();

  nim = NULL;
}

/* ----------------------------------------------------------------------*/
void Brain::destroy() {
  memory->destroy(nv);
  memory->destroy(nvl);
  memory->destroy(npart);

  memory->destroy(boxlo);
  memory->destroy(boxhi);
  memory->destroy(lbox);

  memory->destroy(xlo);
  memory->destroy(xhi);

  memory->destroy(x);
  memory->destroy(type);
  memory->destroy(group);
  memory->destroy(is_loc);
  memory->destroy(tag);

  memory->destroy(num_neigh);
  memory->destroy(neigh);

  memory->destroy(agent);
  memory->destroy(deriv);

  if(nim)
    nifti_image_free(nim);

}
