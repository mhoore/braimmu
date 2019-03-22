#include <mpi.h>
#include <vector>
#include "math.h"

#include "pointers.h"
#include "brain.h"
#include "init.h"

#include <stdint.h>

using namespace std;
using namespace brain_NS;

/* ---------------------------------------------------------------------- */
Init::Init() {
}

/* ----------------------------------------------------------------------*/
Init::~Init() {
  int i;

  for (i=0; i<mri_arg.size(); i++) {
    mri_arg[i].clear();
  }

  mri_arg.clear();
}

/* ----------------------------------------------------------------------*/
void Init::setup(Brain *brn) {
  int me = brn->me;

  read_mri(brn);

  if (!me)
    printf("Setup: setting voxels ... \n");

  boundaries(brn);
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

  set_parameters(brn);

  if (!me)
    printf("Setting communications ... \n");
  brn->comm->comm_init(brn);

  if (!me)
    printf("Setting neighbors ... \n");
  neighbor(brn);

}

/* ----------------------------------------------------------------------
 * Read all mri files, respectively in the order of input
 * ----------------------------------------------------------------------*/
void Init::read_mri(Brain *brn) {
  int i,j;
  int narg;

  // set the main nifti image based on the first mri input
  brn->nim = nifti_image_read(mri_arg[0][1].c_str(),1);

  for (i=0; i<mri_arg.size(); i++) {
    narg = mri_arg[i].size();

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
             mri_arg[i][0].compare("csf")) {
      printf("Error: read_mri unknown keyword. keywords: wm, gm, csf, all, restart. \n");
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

    if (!brn->me) {
      printf("##################### \n");
      printf("NIFTI image %s is read. \n", mri_arg[i][1].c_str());
      printf("NIFTI image properties: ");
      printf("ndim = %i \n", nim_tmp->ndim);
      for (int i=1; i<8; i++)
        printf("dim[%i] = %i, pixdim[i] = %g \n",
               i,nim_tmp->dim[i],i,nim_tmp->pixdim[i]);
      printf("nvox = %lli \n", nim_tmp->nvox);
      printf("nbyper = %i \n", nim_tmp->nbyper);
      printf("datatype = %i \n", nim_tmp->datatype);

      printf("calmin = %g, calmax = %g \n", nim_tmp->cal_min, nim_tmp->cal_max);
      printf("toffset = %g \n", nim_tmp->toffset);

      printf("xyz_units = %i, time_units = %i \n", nim_tmp->xyz_units, nim_tmp->time_units);
      printf("nifti_type = %i \n", nim_tmp->nifti_type);

      printf("intent_code = %i \n", nim_tmp->intent_code);

      printf("description: %s \n", nim_tmp->descrip);
      printf("##################### \n");
    }

    nifti_image_free(nim_tmp);
  }

}

/* ----------------------------------------------------------------------
 * Set boundaries, find the total number of voxels
 * ----------------------------------------------------------------------*/
void Init::boundaries(Brain *brn) {
  int i;

  double vlen = brn->vlen;

  int *npart = brn->npart;
  int *nv = brn->nv;

  double *boxlo = brn->boxlo;
  double *boxhi = brn->boxhi;
  double *lbox = brn->lbox;

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
    for (i=0; i<3; i++) {
      lbox[i] = boxhi[i] - boxlo[i];
      nv[i] = static_cast<int>(lbox[i] / vlen);
    }
  }

  brn->nvoxel = nv[0] * nv[1] * nv[2];

}

/* ----------------------------------------------------------------------
 * Set local and ghost voxels for each partition
 * ----------------------------------------------------------------------*/
void Init::voxels(Brain *brn, int allocated) {
  int i,j,k;

  tagint nvoxel;
  int nlocal,nghost,nall;
  double pos[3];

  double vlen = brn->vlen;

  int *nv = brn->nv;

  double *boxlo = brn->boxlo;

  double *xlo = brn->xlo;
  double *xhi = brn->xhi;

  // find nlocal and nghost
  nlocal = nghost = 0;
  for (i=0; i<nv[0]; i++) {
    pos[0] = boxlo[0] + (0.5 + i) * vlen;
    for (j=0; j<nv[1]; j++) {
      pos[1] = boxlo[1] + (0.5 + j) * vlen;
      for (k=0; k<nv[2]; k++) {
        pos[2] = boxlo[2] + (0.5 + k) * vlen;

        if (pos[0] >= xlo[0] && pos[0] < xhi[0]
         && pos[1] >= xlo[1] && pos[1] < xhi[1]
         && pos[2] >= xlo[2] && pos[2] < xhi[2])
          nlocal++;
        //else if (pos[0] >= xlo[0] - vlen && pos[0] < xhi[0] + vlen
        //      && pos[1] >= xlo[1] - vlen && pos[1] < xhi[1] + vlen
        //      && pos[2] >= xlo[2] - vlen && pos[2] < xhi[2] + vlen)
        else if ( (pos[0] >= xlo[0] - vlen && pos[0] < xhi[0] + vlen
                && pos[1] >= xlo[1] && pos[1] < xhi[1]
                && pos[2] >= xlo[2] && pos[2] < xhi[2])
               || (pos[1] >= xlo[1] - vlen && pos[1] < xhi[1] + vlen
                && pos[0] >= xlo[0] && pos[0] < xhi[0]
                && pos[2] >= xlo[2] && pos[2] < xhi[2])
               || (pos[2] >= xlo[2] - vlen && pos[2] < xhi[2] + vlen
                && pos[0] >= xlo[0] && pos[0] < xhi[0]
                && pos[1] >= xlo[1] && pos[1] < xhi[1]) )
          nghost++;
      }
    }
  }

  nall = nlocal + nghost;
  brn->nlocal = nlocal;
  brn->nghost = nghost;
  brn->nall = nall;

  /// allocations
  allocations(brn, allocated);

  double **x = brn->x;
  tagint *tag = brn->tag;
  int *map = brn->map;

  for (i=0; i<brn->nvoxel; i++)
    map[i] = -1;

  // setup voxel positions, tags, and mapping from tag to id
  nlocal = nvoxel = 0;
  nghost = brn->nlocal;
  for (i=0; i<nv[0]; i++) {
    pos[0] = boxlo[0] + (0.5 + i) * vlen;
    for (j=0; j<nv[1]; j++) {
      pos[1] = boxlo[1] + (0.5 + j) * vlen;
      for (k=0; k<nv[2]; k++) {
        pos[2] = boxlo[2] + (0.5 + k) * vlen;

        if (pos[0] >= xlo[0] && pos[0] < xhi[0]
         && pos[1] >= xlo[1] && pos[1] < xhi[1]
         && pos[2] >= xlo[2] && pos[2] < xhi[2]) {
          x[nlocal][0] = pos[0];
          x[nlocal][1] = pos[1];
          x[nlocal][2] = pos[2];

          tag[nlocal] = nvoxel;
          map[tag[nlocal]] = nlocal;
          nlocal++;
        }
        //else if (pos[0] >= xlo[0] - vlen && pos[0] < xhi[0] + vlen
        //      && pos[1] >= xlo[1] - vlen && pos[1] < xhi[1] + vlen
        //      && pos[2] >= xlo[2] - vlen && pos[2] < xhi[2] + vlen) {
        else if ( (pos[0] >= xlo[0] - vlen && pos[0] < xhi[0] + vlen
                && pos[1] >= xlo[1] && pos[1] < xhi[1]
                && pos[2] >= xlo[2] && pos[2] < xhi[2])
               || (pos[1] >= xlo[1] - vlen && pos[1] < xhi[1] + vlen
                && pos[0] >= xlo[0] && pos[0] < xhi[0]
                && pos[2] >= xlo[2] && pos[2] < xhi[2])
               || (pos[2] >= xlo[2] - vlen && pos[2] < xhi[2] + vlen
                && pos[0] >= xlo[0] && pos[0] < xhi[0]
                && pos[1] >= xlo[1] && pos[1] < xhi[1]) ) {

          x[nghost][0] = pos[0];
          x[nghost][1] = pos[1];
          x[nghost][2] = pos[2];

          tag[nghost] = nvoxel;
          map[tag[nghost]] = nghost;
          nghost++;
        }

        nvoxel++;
      }
    }
  }

  // set topology based on mri input if it exists.
  mri_topology(brn,brn->nim);

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
  //if (brn->nlocal < nvl[0]*nvl[1]*nvl[2]) {
  //  printf("Error333: voxels setup went wrong. proc %i: nlocal = %i, nvl = %i %i %i \n",
  //         brn->me,brn->nlocal,nvl[0],nvl[1],nvl[2]);
  //  exit(1);
  //}
  //printf("Proc %i: %li %li Here2 \n",brn->me, nlocal, nlocal*dsize);
  ////////////////

}

/* ----------------------------------------------------------------------
 * Set neighbors for each voxel: each voxel has 6 neighbors for interaction.
 * The information of only 3 neighbors are saved in each voxel, the other
 * three neighbors are automatically taken into account from other voxels.
 * ----------------------------------------------------------------------*/
void Init::neighbor(Brain *brn) {
  tagint ii,jj,kk,itag;
  int i,j;
  //double dr,delx,dely,delz;

  int nlocal = brn->nlocal;
  //int nall = brn->nall;

  double **x = brn->x;
  int *map = brn->map;

  double vlen_1 = brn->vlen_1;
  double *boxlo = brn->boxlo;

  brn->num_neigh_max = 3;

  create(brn->num_neigh,nlocal,"num_neigh");
  create(brn->neigh,nlocal,brn->num_neigh_max,"neigh");

  int *num_neigh = brn->num_neigh;
  int **neigh = brn->neigh;

  // only one voxel has the info of the neighboring
  for (i=0; i<nlocal; i++) {
    num_neigh[i] = 0;

    ii = (int) ( (x[i][0] - boxlo[0]) * vlen_1 - 0.5 );
    jj = (int) ( (x[i][1] - boxlo[1]) * vlen_1 - 0.5 );
    kk = (int) ( (x[i][2] - boxlo[2]) * vlen_1 - 0.5 );

    itag = find_me(brn,ii+1,jj,kk);
    if (itag != -1) {
      j = map[itag];
      neigh[i][num_neigh[i]] = j;
      num_neigh[i]++;
    }

    itag = find_me(brn,ii,jj+1,kk);
    if (itag != -1) {
      j = map[itag];
      neigh[i][num_neigh[i]] = j;
      num_neigh[i]++;
    }

    itag = find_me(brn,ii,jj,kk+1);
    if (itag != -1) {
      j = map[itag];
      neigh[i][num_neigh[i]] = j;
      num_neigh[i]++;
    }

    /// DEBUG
    ////////////////
    if (brn->tag[i] != find_me(brn,ii,jj,kk)) {
      printf("Error: neighboring. proc %i: " TAGINT_FORMAT " " TAGINT_FORMAT " \n",
             brn->me,brn->tag[i],find_me(brn,ii,jj,kk));
      exit(1);
    }
   ///////
    //for (j=0; j<nall; j++) {
      //delx = (x[j][0] - x[i][0]) * vlen_1;
      //if (delx < 0.0) continue;
      //if (delx > 1.1) continue;

      //dely = (x[j][1] - x[i][1]) * vlen_1;
      //if (dely < 0.0) continue;
      //if (dely > 1.1) continue;

      //delz = (x[j][2] - x[i][2]) * vlen_1;
      //if (delz < 0.0) continue;
      //if (delz > 1.1) continue;

      //dr = sqrt(delx*delx + dely*dely + delz*delz);
      //if (dr > 1.1) continue;

      //neigh[i][num_neigh[i]] = j;
      //num_neigh[i]++;
      //}
  }

}

/* ----------------------------------------------------------------------
 * Find the tag of a voxel from its coordinates i,j,k
 * ----------------------------------------------------------------------*/
tagint Init::find_me(Brain *brn, int i, int j, int k) {
  int *nv = brn->nv;
  tagint itag;

  if (i < 0 || i >= nv[0])
    return -1;
  if (j < 0 || j >= nv[1])
    return -1;
  if (k < 0 || k >= nv[2])
    return -1;

  itag = k + nv[2] * (j + nv[1]*i);

  return itag;

}

/* ----------------------------------------------------------------------*/
void Init::allocations(Brain *brn, int allocated) {

  if (allocated) {
    destroy(brn->type);
    destroy(brn->group);

    destroy(brn->x);
    destroy(brn->tag);

    destroy(brn->agent);
    //destroy(brn->grad);
  }

  create(brn->type,brn->nall,"type");
  create(brn->group,brn->nall,"group");

  create(brn->x,brn->nall,3,"x");
  create(brn->tag,brn->nall,"tag");

  if (!allocated)
    create(brn->map,brn->nvoxel,"map");

  create(brn->agent,num_agents,brn->nall,2,"agent");
  //create(brn->grad,num_agents,brn->nall,3,"grad");

  // set initial values
  for (int ag_id=0; ag_id<num_agents; ag_id++) {
    if (brn->init_val[ag_id] >= 0.0)
      for (int i=0; i<brn->nall; i++)
        brn->agent[ag_id][i][0] = brn->init_val[ag_id];
  }

}

/* ----------------------------------------------------------------------
 * Set constant global simulation parameters.
 * ----------------------------------------------------------------------*/
void Init::set_parameters(Brain *brn) {
  brn->D_sAb = brn->diff_sAb * brn->vlen_2;
  brn->D_mic = brn->diff_mic * brn->vlen_2;
  brn->cs = brn->sens_s * brn->vlen_2;
  brn->cf = brn->sens_f * brn->vlen_2;
}

/* ----------------------------------------------------------------------
 * Define the boundaries of the system based on the mri nifti image (.nii)
 * ----------------------------------------------------------------------*/
int Init::mri_boundaries(Brain *brn, nifti_image *nim) {
  if (!nim)
    return 0;

  int i;

  double vlen = brn->vlen;

  int *nv = brn->nv;

  double *boxlo = brn->boxlo;
  double *boxhi = brn->boxhi;
  double *lbox = brn->lbox;

  double conver_fac = 1.0;
  if (nim->xyz_units == NIFTI_UNITS_METER)
    conver_fac = 1.e6;
  else if (nim->xyz_units == NIFTI_UNITS_MM)
    conver_fac = 1.e3;
  else if (nim->xyz_units == NIFTI_UNITS_MICRON)
    conver_fac = 1.0;

  for (i=0; i<3; i++) {
    lbox[i] = nim->dim[i+1] * nim->pixdim[i+1];
    lbox[i] *= conver_fac;
    boxlo[i] = -0.5 * lbox[i];
    boxhi[i] = 0.5 * lbox[i];

    nv[i] = static_cast<int>(lbox[i] / vlen);

  }

  if (!brn->me) {
    printf("Number of voxels in the original mri file: nx ny nz \n");
    printf("%i %i %i \n", nim->dim[1], nim->dim[2], nim->dim[3]);
    printf("Number of voxels assigned: nx ny nz \n");
    printf("%i %i %i \n", nv[0], nv[1], nv[2]);
  }

  return 1;

}

/* ----------------------------------------------------------------------
 * Define the system topology based on the mri nifti image (.nii)
 * ----------------------------------------------------------------------*/
void Init::mri_topology(Brain *brn, nifti_image *nim) {
  if (!nim)
    return;

  int i,j,k,h;

  int nlocal = brn->nlocal;
  int nall = brn->nall;

  int *map = brn->map;

  int *type = brn->type;
  double ***agent = brn->agent;

  double vlen_1 = brn->vlen_1;

  int ii,jj,kk,vid;
  double dum;
  tagint itag;

  double conver_fac = 1.0;
  if (nim->xyz_units == NIFTI_UNITS_METER)
    conver_fac = 1.e6;
  else if (nim->xyz_units == NIFTI_UNITS_MM)
    conver_fac = 1.e3;
  else if (nim->xyz_units == NIFTI_UNITS_MICRON)
    conver_fac = 1.0;

  // set all voxel types as EMP_type
  for (i=0; i<nall; i++)
    type[i] = EMP_type;

  /* -------------------------------------------------------
   * set from restart
   * ------------------------------------------------------- */
  if (!mri_arg[0][0].compare("restart")) {
    double **v_prop;
    int **n_prop;

    brn->memory->create(v_prop,nall,nim->dim[5],"init:v_prop");
    brn->memory->create(n_prop,nall,nim->dim[5],"init:n_prop");

    for (i=0; i<nall; i++)
      for (j=0; j<nim->dim[5]; j++) {
        v_prop[i][j] = 0.0;
        n_prop[i][j] = 0;
      }

    // set pointers to data
    if (nim->datatype == DT_UINT8)
      ptr8 = (uint8_t *) nim->data;
    else if (nim->datatype == DT_INT16)
      ptr16 = (int16_t *) nim->data;
    else if (nim->datatype == DT_INT32)
      ptr32 = (int32_t *) nim->data;
    else if (nim->datatype == DT_FLOAT32)
      ptrf = (float *) nim->data;
    else if (nim->datatype == DT_FLOAT64)
      ptrd = (double *) nim->data;
    else {
      printf("Error: nifti file data type cannot be read. datatype=%i . \n", nim->datatype);
      exit(1);
    }

    int c = 0;
    for (h=0; h<nim->dim[5]; h++) {
      for (k=0; k<nim->dim[3]; k++) {
        dum = nim->pixdim[3] * k * conver_fac;
        kk = (int) (dum * vlen_1);

        for (j=0; j<nim->dim[2]; j++) {
          dum = nim->pixdim[2] * j * conver_fac;
          jj = (int) (dum * vlen_1);

          for (i=0; i<nim->dim[1]; i++) {
            dum = nim->pixdim[1] * i * conver_fac;
            ii = (int) (dum * vlen_1);

            itag = find_me(brn,ii,jj,kk);
            if (itag == -1) {
              //printf("Warning: a tag cannot be assigned for the voxels of the mri file. \n");
              c++;
              continue;
            }
            vid = map[itag];

            // if it is in the partition
            if (vid != -1) {
              if (nim->datatype == DT_UINT8)
                v_prop[vid][h] += (double) ptr8[c];
              else if (nim->datatype == DT_INT16)
                v_prop[vid][h] += (double) ptr16[c];
              else if (nim->datatype == DT_INT32)
                v_prop[vid][h] += (double) ptr32[c];
              else if (nim->datatype == DT_FLOAT32)
                v_prop[vid][h] += (double) ptrf[c];
              else if (nim->datatype == DT_FLOAT64)
                v_prop[vid][h] += (double) ptrd[c];

              n_prop[vid][h]++;
            }
            c++;
          }
        }
      }
    }

    string str = nim->descrip;
    vector<string> arg;
    arg.clear();

    istringstream buf(str);

    for(string word; buf >> word;)
      arg.push_back(word);

    int narg = arg.size();

    if (narg != nim->dim[5]){
      printf("Error: mri file contents do not match its description. \n");
      exit(1);
    }

    // set voxel properties based on the nifti_image data
    for (i=0; i<nlocal; i++) {
      for (j=0; j<nim->dim[5]; j++) {

        if (n_prop[i][j] > 0)
          v_prop[i][j] /= n_prop[i][j];

        if (!arg[j].compare("type"))
          type[i] = (int) round(v_prop[i][j]);
        else if (brn->input->find_agent(arg[j]) >= 0)
          agent[brn->input->find_agent(arg[j])][i][0] = v_prop[i][j];
        else {
          printf("Error: mri file content cannot be assigned. arg = %s \n", arg[j].c_str());
          exit(1);
        }
      }
    }

    brn->memory->destroy(v_prop);
    brn->memory->destroy(n_prop);

  }

  /* -------------------------------------------------------
   * if there is no restart, go through the nifti_image data
   * and find the corresponding voxel.
   * ------------------------------------------------------- */
  else {
    if (nim->ndim != 3) {
      printf("Error: nifti file should have 3 dimensions. ndim = %i \n",nim->ndim);
      exit(1);
    }

    double max_val, thres_val;

    double *v_prop;
    int *n_prop;

    brn->memory->create(v_prop,nall,"init:v_prop");
    brn->memory->create(n_prop,nall,"init:n_prop");

    // go though all mri files
    for (int tis=0; tis<mri_arg.size(); tis++) {
      nifti_image *nim_tmp = NULL;

      nim_tmp = nifti_image_read(mri_arg[tis][1].c_str(),1);
      thres_val = stof(mri_arg[tis][2]);
      max_val = stof(mri_arg[tis][3]);

      // set pointers to data
      if (nim_tmp->datatype == DT_UINT8)
        ptr8 = (uint8_t *) nim_tmp->data;
      else if (nim_tmp->datatype == DT_INT16)
        ptr16 = (int16_t *) nim_tmp->data;
      else if (nim_tmp->datatype == DT_INT32)
        ptr32 = (int32_t *) nim_tmp->data;
      else if (nim_tmp->datatype == DT_FLOAT32)
        ptrf = (float *) nim_tmp->data;
      else if (nim_tmp->datatype == DT_FLOAT64)
        ptrd = (double *) nim_tmp->data;
      else {
        printf("Error: nifti file data type cannot be read. datatype=%i . \n", nim_tmp->datatype);
        exit(1);
      }

      for (i=0; i<nall; i++) {
        v_prop[i] = 0.0;
        n_prop[i] = 0;
      }

      int c = 0;
      for (k=0; k<nim_tmp->dim[3]; k++) {
        dum = nim_tmp->pixdim[3] * k * conver_fac;
        kk = (int) (dum * vlen_1);

        for (j=0; j<nim_tmp->dim[2]; j++) {
          dum = nim_tmp->pixdim[2] * j * conver_fac;
          jj = (int) (dum * vlen_1);

          for (i=0; i<nim_tmp->dim[1]; i++) {
            dum = nim_tmp->pixdim[1] * i * conver_fac;
            ii = (int) (dum * vlen_1);

            itag = find_me(brn,ii,jj,kk);
            if (itag == -1) {
              //printf("Warning: a tag cannot be assigned for the voxels of the mri file. \n");
              c++;
              continue;
            }
            vid = map[itag];

            double uint8coef = 1.0 / ((double) UINT8_MAX);
            double int16coef = 1.0 / ((double) INT16_MAX - (double) INT16_MIN);
            double int32coef = 1.0 / ((double) INT32_MAX - (double) INT32_MIN);
            // if it is in the partition
            if (vid != -1) {
              if (nim_tmp->datatype == DT_UINT8)
                v_prop[vid] += (double) ptr8[c] * uint8coef;
              else if (nim_tmp->datatype == DT_INT16)
                v_prop[vid] += ((double) (ptr16[c] - INT16_MIN)) * int16coef;
              else if (nim_tmp->datatype == DT_INT32)
                v_prop[vid] += ((double) (ptr32[c] - INT32_MIN)) * int32coef;
              else if (nim_tmp->datatype == DT_FLOAT32)
                v_prop[vid] += (double) ptrf[c];
              else if (nim_tmp->datatype == DT_FLOAT64)
                v_prop[vid] += (double) ptrd[c];
              n_prop[vid]++;

              //printf("HERE: proc %i: vid = %i, v_prop[vid] = %g, ptr16[c] = %i, c = %i, %i \n",
                     //brn->me, vid, v_prop[vid], ptr16[c], c, nim_tmp->datatype);

            }

            c++;
          }
        }
      }

      // set voxel properties based on the nifti_image data
      for (i=0; i<nlocal; i++) {
        if (n_prop[i] > 0)
          v_prop[i] /= n_prop[i];

        double coef = 1.0 / (1.0 - thres_val);
        // criteria based on mri file
        if (!mri_arg[tis][0].compare("all")) {
          if (v_prop[i] <= 0) {
            type[i] = EMP_type;
            for (int ag_id=0; ag_id<num_agents; ag_id++)
              agent[ag_id][i][0] = 0.0;
          }

          else if (v_prop[i] < thres_val) {
            type[i] = CSF_type;
            for (int ag_id=0; ag_id<num_agents; ag_id++)
              agent[ag_id][i][0] = 0.0;
          }

          else if (v_prop[i] > thres_val) {
            type[i] = GM_type;
            agent[neu][i][0] = (v_prop[i] - thres_val) * coef * max_val;
          }
        }

        else if (!mri_arg[tis][0].compare("wm")) {
          if (v_prop[i] > thres_val) {
            type[i] = WM_type;
            agent[neu][i][0] += (v_prop[i] - thres_val) * coef * max_val;
          }
        }

        else if (!mri_arg[tis][0].compare("gm")) {
          if (v_prop[i] > thres_val) {
            type[i] = GM_type;
            agent[neu][i][0] += (v_prop[i] - thres_val) * coef * max_val;
          }
        }

        else if (!mri_arg[tis][0].compare("csf")) {
          if (v_prop[i] > thres_val) {
            type[i] = CSF_type;
            for (int ag_id=0; ag_id<num_agents; ag_id++)
              agent[ag_id][i][0] = 0.0;
          }
        }

      }

      nifti_image_free(nim_tmp);

    }

    brn->memory->destroy(v_prop);
    brn->memory->destroy(n_prop);

  }

  // set all values to zero for EMT_type voxels
  for (i=0; i<nall; i++)
    if (type[i] == EMP_type) {
      for (int ag_id=0; ag_id<num_agents; ag_id++)
        agent[ag_id][i][0] = 0.0;
    }

}

/* ----------------------------------------------------------------------
 * Using map function rather than map array (pointer) for less memory usage
 * purpose. Otherwise, the map array is a faster process.
 * ----------------------------------------------------------------------*/
/*
int map(Brain*, tagint); // mapping from tag to id

// find approximate nvl
double dx[3];
for (i=0; i<3; i++) {
  dx[i] = xhi[i] - xlo[i];
  nvl[i] = (int) (dx[i] * brn->vlen_1);
  if ((0.5 + nvl[i]) * vlen < dx[i])
    nvl[i] += 1;
}

int Init::map(Brain *brn, tagint itag) {
  tagint dumt;
  int i,j,k,istart;

  double pos[3],dx[3];

  double vlen = brn->vlen;
  double vlen_1 = brn->vlen_1;

  int *nvl = brn->nvl;

  double *boxlo = brn->boxlo;
  double *xlo = brn->xlo;
  double *xhi = brn->xhi;

  int *nv = brn->nv;

  tagint *tag = brn->tag;

  i = (int) (itag / (nv[1]*nv[2]));
  dumt = itag % (nv[1]*nv[2]);
  j = (int) (dumt / nv[2]);
  k = dumt % nv[2];

  pos[0] = boxlo[0] + (0.5 + i) * vlen;
  pos[1] = boxlo[1] + (0.5 + j) * vlen;
  pos[2] = boxlo[2] + (0.5 + k) * vlen;

  if (pos[0] < xhi[0] && pos[1] < xhi[1] && pos[2] < xhi[2]) {
    dx[0] = pos[0] - xlo[0] + 0.5;
    dx[1] = pos[1] - xlo[1] + 0.5;
    dx[2] = pos[2] - xlo[2] + 0.5;

    i = (int) (dx[0] * vlen_1);
    j = (int) (dx[1] * vlen_1);
    k = (int) (dx[2] * vlen_1);

    istart = i*nvl[1]*nvl[2] + j*nvl[2] + k;
  }
  else if (pos[0] < xhi[0] + vlen && pos[1] < xhi[1] + vlen && pos[2] < xhi[2] + vlen)
     istart = brn->nlocal;
  else
    return -1;

  for (i=istart; i<brn->nall; i++)
    if (tag[i] == itag)
      return i;

  return -1;

}
*/
