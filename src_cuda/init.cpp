#include "virtualbrain.h"
#include "init.h"

using namespace std;

/* ---------------------------------------------------------------------- */
Init::Init() {
}

/* ----------------------------------------------------------------------*/
Init::~Init() {
}

/* ----------------------------------------------------------------------*/
void Init::setup(VirtualBrain *brn) {
  int me = brn->me;

  if (mri_arg.size() > 0)
    read_mri(brn);

  if (!me)
    printf("Setup: setting voxels ... \n");

  boundaries(brn);

  brn->set_parameters();

  brn->comm->partition(brn);
  voxels(brn,0);
  brn->region->apply_regions(brn);

  int itr = 0;
  if (!me && brn->comm->b_itr)
    printf("Balancing workload on the partitions ... \n");
  while (itr < brn->comm->b_itr) {
    // balance xlo and xhi in dimension x[i]
    brn->comm->balance(brn);
    voxels(brn,1);
    brn->region->apply_regions(brn);
    itr++;
  }

  if (!me)
    printf("Setting communications ... \n");
  brn->comm->comm_init(brn);

}

/* ----------------------------------------------------------------------
 * Read all mri files, respectively in the order of input
 * ----------------------------------------------------------------------*/
void Init::read_mri(VirtualBrain *brn) {
  // set the main nifti image based on the first mri input
  brn->nim = nifti_image_read(mri_arg[0][1].c_str(),1);

  for (int i=0; i<mri_arg.size(); i++) {
    int narg = mri_arg[i].size();

    if (narg < 3){
      printf("Error: read_mri has fewer arguments than required. \n");
      exit(1);
    }

    if (!mri_arg[i][0].compare("restart")) {
      if (mri_arg.size() > 1) {
        printf("Error: read_mri restart should be used alone. Remove other read_mri commands. \n");
        exit(1);
      }
    }

    else if (!mri_arg[i][0].compare("all")) {
      if (mri_arg.size() > 1) {
        printf("Error: read_mri all should be used alone. Remove other read_mri commands. \n");
        exit(1);
      }
    }

    else if (mri_arg[i][0].compare("wm") &&
             mri_arg[i][0].compare("gm") &&
             mri_arg[i][0].compare("csf") &&
             mri_arg[i][0].compare("group") &&
             mri_arg[i][0].compare("rgb") ) {
      printf("Error: read_mri unknown keyword. keywords: wm, gm, csf, group, rgb, all, restart. \n");
      exit(1);
    }

    nifti_image *nim_tmp;
    nim_tmp = NULL;

    /* read the nifti header, and optionally image data */
    nim_tmp = nifti_image_read(mri_arg[i][1].c_str(),0);
    if(!nim_tmp) {
      printf("Error: read_mri for %s - nifti image %s not readable. \n",
             mri_arg[i][0].c_str(), mri_arg[i][1].c_str());
      exit(1);
    }

    if (!brn->me)
      print_mri_properties(brn,nim_tmp,mri_arg[i][1]);

    nifti_image_free(nim_tmp);
  }

}

/* ----------------------------------------------------------------------
 * Print important properties of MRI data file
 * ----------------------------------------------------------------------*/
void Init::print_mri_properties(VirtualBrain *brn, nifti_image *nim, string fname) {
  ofstream logfile;
  logfile.open (flog, ios::app);

  char buf[256];
  sprintf(buf,"##################### ");
  cout << buf << endl;
  logfile << buf << endl;

  sprintf(buf,"Reading NIFTI image %s: ", fname.c_str());
  cout << buf << endl;
  logfile << buf << endl;

  sprintf(buf,"ndim = %i ", nim->ndim);
  cout << buf << endl;
  logfile << buf << endl;

  for (int i=1; i<=nim->ndim; i++) {
    sprintf(buf,"dim[%i] = %i, pixdim[i] = %g ",
           i,nim->dim[i],i,nim->pixdim[i]);
    cout << buf << endl;
    logfile << buf << endl;
  }

  sprintf(buf,"nvox = %lli ", nim->nvox);
  cout << buf << endl;
  logfile << buf << endl;

  sprintf(buf,"nbyper = %i ", nim->nbyper);
  cout << buf << endl;
  logfile << buf << endl;

  sprintf(buf,"datatype = %i ", nim->datatype);
  cout << buf << endl;
  logfile << buf << endl;

  sprintf(buf,"xyz_units = %i, time_units = %i ", nim->xyz_units, nim->time_units);
  cout << buf << endl;
  logfile << buf << endl;

  sprintf(buf,"nifti_type = %i ", nim->nifti_type);
  cout << buf << endl;
  logfile << buf << endl;

  sprintf(buf,"description: %s ", nim->descrip);
  cout << buf << endl;
  logfile << buf << endl;

  sprintf(buf,"##################### ");
  cout << buf << endl;
  logfile << buf << endl;

  logfile.close();

}

/* ----------------------------------------------------------------------
 * Set boundaries, find the total number of voxels
 * ----------------------------------------------------------------------*/
void Init::boundaries(VirtualBrain *brn) {
  double vlen = brn->vlen;

  auto &npart = brn->npart;
  auto &nv = brn->nv;

  auto &boxlo = brn->boxlo;
  auto &boxhi = brn->boxhi;
  auto &lbox = brn->lbox;

  if (brn->nproc != npart[0] * npart[1] * npart[2]) {
    printf("Error: partitions mismatch number of processors. \n");
    exit(1);
  }

  brn->vlen_1 = 1.0 / vlen;
  brn->vlen_2 = brn->vlen_1 * brn->vlen_1;
  brn->vvol = vlen * vlen * vlen * 1.0e-15; // voxel volume in mL
  brn->vvol_1 = 1.0 / brn->vvol; // 1/mL

  /* Set boundaries based on input mri file if input mri exists. */
  if (!mri_boundaries(brn,brn->nim)) {
    for (int i=0; i<3; i++) {
      lbox[i] = boxhi[i] - boxlo[i];
      nv[i] = static_cast<int>( round(lbox[i] / vlen) );
    }
  }

  brn->nvoxel = nv[0] * nv[1] * nv[2];

}

/* ----------------------------------------------------------------------
 * Set local and ghost voxels for each partition
 * ----------------------------------------------------------------------*/
void Init::voxels(VirtualBrain *brn, bool allocated) {
  tagint nvoxel;
  int nlocal,nghost,nall;
  double pos[3];

  double vlen = brn->vlen;

  auto &nv = brn->nv;

  auto &boxlo = brn->boxlo;

  auto &xlo = brn->xlo;
  auto &xhi = brn->xhi;

  // find nlocal and nghost
  nlocal = nghost = 0;
  for (int k=-1; k<nv[2]+1; k++) {
    pos[2] = boxlo[2] + (0.5 + k) * vlen;
    for (int j=-1; j<nv[1]+1; j++) {
      pos[1] = boxlo[1] + (0.5 + j) * vlen;
      for (int i=-1; i<nv[0]+1; i++) {
        pos[0] = boxlo[0] + (0.5 + i) * vlen;

        if ( pos[0] >= xlo[0] && pos[0] < xhi[0]
          && pos[1] >= xlo[1] && pos[1] < xhi[1]
          && pos[2] >= xlo[2] && pos[2] < xhi[2] )
          nlocal++;
        else if ( pos[0] >= xlo[0] - vlen && pos[0] < xhi[0] + vlen
               && pos[1] >= xlo[1] - vlen && pos[1] < xhi[1] + vlen
               && pos[2] >= xlo[2] - vlen && pos[2] < xhi[2] + vlen )
          nghost++;
      }
    }
  }

  nall = nlocal + nghost;
  brn->nlocal = nlocal;
  brn->nghost = nghost;
  brn->nall = nall;

  if (!allocated) {
    maparr.clear();
    maparr.resize(brn->nvoxel);
  }

  brn->allocations();

  auto &x = brn->x;
  auto &tag = brn->tag;

  auto &is_loc = brn->is_loc;

  for (tagint i=0; i<brn->nvoxel; i++)
    maparr[i] = -1;

  // setup voxel positions, tags, and mapping from tag to id
  nvoxel = nall = 0;
  for (int k=-1; k<nv[2]+1; k++) {
    pos[2] = boxlo[2] + (0.5 + k) * vlen;
    for (int j=-1; j<nv[1]+1; j++) {
      pos[1] = boxlo[1] + (0.5 + j) * vlen;
      for (int i=-1; i<nv[0]+1; i++) {
        pos[0] = boxlo[0] + (0.5 + i) * vlen;

        if ( pos[0] >= xlo[0] && pos[0] < xhi[0]
          && pos[1] >= xlo[1] && pos[1] < xhi[1]
          && pos[2] >= xlo[2] && pos[2] < xhi[2] ) {
          x[0][nall] = pos[0];
          x[1][nall] = pos[1];
          x[2][nall] = pos[2];

          tag[nall] = nvoxel;
          is_loc[nall] = 1;
          maparr[nvoxel] = nall;
          nall++;
        }

        else if ( pos[0] >= xlo[0] - vlen && pos[0] < xhi[0] + vlen
               && pos[1] >= xlo[1] - vlen && pos[1] < xhi[1] + vlen
               && pos[2] >= xlo[2] - vlen && pos[2] < xhi[2] + vlen ) {
          x[0][nall] = pos[0];
          x[1][nall] = pos[1];
          x[2][nall] = pos[2];

          if (k>=0 && k<nv[2]
           && j>=0 && j<nv[1]
           && i>=0 && i<nv[0]) {
            tag[nall] = nvoxel;
            maparr[nvoxel] = nall;
          }
          else
            tag[nall] = -1;

          is_loc[nall] = 0;
          nall++;
        }

        if (k>=0 && k<nv[2]
         && j>=0 && j<nv[1]
         && i>=0 && i<nv[0])
          nvoxel++;
      }
    }
  }

  // set nvl
  for (int i=0; i<ndim; i++)
    brn->nvl[i] =
        static_cast<int>( round( (x[i][brn->nall - 1] - x[i][0]) * brn->vlen_1 ) ) - 1;

  // set topology based on mri input if it exists.
  brn->mri_topology(brn->nim);

  /// DEBUG
  ////////////////
  if (nvoxel != brn->nvoxel) {
    printf("Error111: voxels setup went wrong. proc %i: nvoxel = %li, brn->nvoxel = %li \n",
           brn->me,nvoxel,brn->nvoxel);
    exit(1);
  }
  tagint nlocal_tmp = nlocal;
  MPI_Allreduce(&nlocal_tmp,&nvoxel,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  if (nvoxel != brn->nvoxel) {
    printf("Error222: voxels setup went wrong. proc %i: nvoxel = %li, brn->nvoxel = %li \n",
          brn->me,nvoxel,brn->nvoxel);
    exit(1);
  }
  ////////////////

}

/* ----------------------------------------------------------------------
 * Find the tag of a voxel from its global coordinates i,j,k
 * ----------------------------------------------------------------------*/
tagint Init::find_tag(VirtualBrain *brn, int i, int j, int k) {
  auto &nv = brn->nv;

  if (i < 0 || i >= nv[0])
    return -1;
  if (j < 0 || j >= nv[1])
    return -1;
  if (k < 0 || k >= nv[2])
    return -1;

  return i + nv[0] * (j + nv[1]*k);

}

/* ----------------------------------------------------------------------
 * Define the boundaries of the system based on the mri nifti image (.nii)
 * ----------------------------------------------------------------------*/
int Init::mri_boundaries(VirtualBrain *brn, nifti_image *nim) {
  if (!nim)
    return 0;

  double vlen = brn->vlen;

  auto &nv = brn->nv;

  auto &boxlo = brn->boxlo;
  auto &boxhi = brn->boxhi;
  auto &lbox = brn->lbox;

  double conver_fac = 1.0;
  if (nim->xyz_units == NIFTI_UNITS_METER)
    conver_fac = 1.e6;
  else if (nim->xyz_units == NIFTI_UNITS_MM)
    conver_fac = 1.e3;
  else if (nim->xyz_units == NIFTI_UNITS_MICRON)
    conver_fac = 1.0;

  for (int i=0; i<3; i++) {
    lbox[i] = nim->dim[i+1] * nim->pixdim[i+1];
    lbox[i] *= conver_fac;
    boxlo[i] = -0.5 * lbox[i];
    boxhi[i] = 0.5 * lbox[i];

    nv[i] = static_cast<int>( round(lbox[i] / vlen) );

  }

  if (!brn->me) {
    printf("Number of voxels in the original mri file: nx ny nz \n");
    printf("%i %i %i \n", nim->dim[1], nim->dim[2], nim->dim[3]);
    printf("Number of voxels assigned: nx ny nz \n");
    printf("%i %i %i \n", nv[0], nv[1], nv[2]);
  }

  return 1;

}
