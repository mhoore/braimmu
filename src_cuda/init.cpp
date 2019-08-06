#include <mpi.h>
#include <vector>
#include "math.h"

#include "pointers.h"
#include "brain.h"
#include "init.h"

#include <stdint.h>

using namespace std;
using namespace brain_NS;

#define PI 3.141592653589793

/* ---------------------------------------------------------------------- */
Init::Init() {
}

/* ----------------------------------------------------------------------*/
Init::~Init() {
}

/* ----------------------------------------------------------------------*/
void Init::setup(Brain *brn) {
  int me = brn->me;

  if (mri_arg.size() > 0)
    read_mri(brn);

  if (!me)
    printf("Setup: setting voxels ... \n");

  boundaries(brn);

  set_parameters(brn);

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
void Init::read_mri(Brain *brn) {
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
void Init::print_mri_properties(Brain *brn, nifti_image *nim, string fname) {
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
void Init::boundaries(Brain *brn) {
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
void Init::voxels(Brain *brn, int allocated) {
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

  /// allocations
  allocations(brn, allocated);

  auto &x = brn->x;
  auto &tag = brn->tag;

  auto &is_loc = brn->is_loc;

  for (tagint i=0; i<brn->nvoxel; i++)
    map[i] = -1;

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
          map[nvoxel] = nall;
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
            map[nvoxel] = nall;
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
  ////////////////

}

/* ----------------------------------------------------------------------
 * Find the tag of a voxel from its global coordinates i,j,k
 * ----------------------------------------------------------------------*/
tagint Init::find_tag(Brain *brn, int i, int j, int k) {
  auto &nv = brn->nv;

  if (i < 0 || i >= nv[0])
    return -1;
  if (j < 0 || j >= nv[1])
    return -1;
  if (k < 0 || k >= nv[2])
    return -1;

  return i + nv[0] * (j + nv[1]*k);

}

/* ----------------------------------------------------------------------*/
void Init::allocations(Brain *brn, int allocated) {

  for (auto &a: brn->x) {
    a.clear();
    a.resize(brn->nall);
  }

  //brn->x.clear();
  //brn->x.resize(brn->nall);

  brn->tag.clear();
  brn->tag.resize(brn->nall);

  brn->type.clear();
  brn->type.resize(brn->nall);

  brn->group.clear();
  brn->group.resize(brn->nall);

  brn->is_loc.clear();
  brn->is_loc.resize(brn->nall);

  for (auto &a: brn->agent) {
    a.clear();
    a.resize(brn->nall);
  }

  for (auto &a: brn->deriv) {
    a.clear();
    a.resize(brn->nall);
  }

  for (auto &a: brn->Dtau) {
    a.clear();
    a.resize(brn->nall);
  }

  if (!allocated) {
    map.clear();
    map.resize(brn->nvoxel);
  }

  // set initial values
  for (int ag_id=0; ag_id<num_agents; ag_id++) {
    if (brn->init_val[ag_id] >= 0.0)
      for (int i=0; i<brn->nall; i++)
        brn->agent[ag_id][i] = brn->init_val[ag_id];
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
  brn->omega_cir = 2.0 * PI / brn->tau_cir;
  brn->Dtau_max = brn->diff_tau * brn->vlen_2;
}

/* ----------------------------------------------------------------------
 * Define the boundaries of the system based on the mri nifti image (.nii)
 * ----------------------------------------------------------------------*/
int Init::mri_boundaries(Brain *brn, nifti_image *nim) {
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

/* ----------------------------------------------------------------------
 * Define the system topology based on the mri nifti image (.nii)
 * ----------------------------------------------------------------------*/
void Init::mri_topology(Brain *brn, nifti_image *nim) {
  if (!nim)
    return;

  int nall = brn->nall;

  auto &tissue = brn->tissue;
  auto &type = brn->type;
  auto &group = brn->group;
  auto &agent = brn->agent;

  auto &Dtau = brn->Dtau;

  double vlen_1 = brn->vlen_1;

  double conver_fac = 1.0;
  if (nim->xyz_units == NIFTI_UNITS_METER)
    conver_fac = 1.e6;
  else if (nim->xyz_units == NIFTI_UNITS_MM)
    conver_fac = 1.e3;
  else if (nim->xyz_units == NIFTI_UNITS_MICRON)
    conver_fac = 1.0;

  // set all voxel types as EMP and groups as 0
  fill(type.begin(),type.end(),tissue[EMP]);
  fill(group.begin(),group.end(),0);
  for (auto &a: Dtau)
    fill(a.begin(),a.end(),brn->Dtau_max);

  /* -------------------------------------------------------
   * set from restart
   * ------------------------------------------------------- */
  if (!mri_arg[0][0].compare("restart")) {
    vector<vector<double>> v_prop(nim->dim[5]);
    for (auto &a: v_prop) {
      a.clear();
      a.resize(nall);
      fill(a.begin(), a.end(), 0.);
    }

    vector<vector<int>> n_prop(nim->dim[5]);
    for (auto &a: n_prop) {
      a.clear();
      a.resize(nall);
      fill(a.begin(), a.end(), 0);
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
    for (int h=0; h<nim->dim[5]; h++) {
      for (int k=0; k<nim->dim[3]; k++) {
        int kk = static_cast<int>( round(nim->pixdim[3] * k * conver_fac * vlen_1) );

        for (int j=0; j<nim->dim[2]; j++) {
          int jj = static_cast<int>( round(nim->pixdim[2] * j * conver_fac * vlen_1) );

          for (int i=0; i<nim->dim[1]; i++) {
            int ii = static_cast<int>( round(nim->pixdim[1] * i * conver_fac * vlen_1) );

            tagint itag = find_tag(brn,ii,jj,kk);
            if (itag == -1) {
              //printf("Warning: a tag cannot be assigned for the voxels of the mri file. \n");
              c++;
              continue;
            }

            int vid = map[itag];

            // if it is in the partition
            if (vid != -1) {
              if (nim->datatype == DT_UINT8)
                v_prop[h][vid] += static_cast<double>(ptr8[c]);
              else if (nim->datatype == DT_INT16)
                v_prop[h][vid] += static_cast<double>(ptr16[c]);
              else if (nim->datatype == DT_INT32)
                v_prop[h][vid] += static_cast<double>(ptr32[c]);
              else if (nim->datatype == DT_FLOAT32)
                v_prop[h][vid] += static_cast<double>(ptrf[c]);
              else if (nim->datatype == DT_FLOAT64)
                v_prop[h][vid] += static_cast<double>(ptrd[c]);

              n_prop[h][vid]++;
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
    for (int i=0; i<nim->dim[5]; i++) {
      for (int j=0; j<nall; j++) {

        if (n_prop[i][j] > 0)
          v_prop[i][j] /= n_prop[i][j];

        if (!arg[i].compare("type"))
          type[j] = static_cast<int>( round(v_prop[i][j]) ); // NOTE: this may lead to errors
        else if (!arg[i].compare("group")) {
          double fractpart, intpart;
          fractpart = modf (v_prop[i][j], &intpart);
          if (fractpart != 0) continue;
          group[j] = static_cast<int>( v_prop[i][j] );
        }
        else if (brn->input->find_agent(arg[i]) >= 0)
          agent[brn->input->find_agent(arg[i])][j] = v_prop[i][j];
        else {
          printf("Error: mri file content cannot be assigned. arg = %s \n", arg[i].c_str());
          exit(1);
        }
      }
    }

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

    vector<double> v_prop(nall);
    vector<int> n_prop(nall);

    vector<vector<double>> rgb_prop(3);
    for (auto &a: rgb_prop) {
      a.clear();
      a.resize(nall);
    }

    // go through all mri files
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
      else if (nim_tmp->datatype == DT_RGB24)
        ptr_rgb = (uint8_t *) nim_tmp->data;
      else {
        printf("Error: nifti file data type cannot be read. datatype=%i . \n", nim_tmp->datatype);
        exit(1);
      }

      fill(v_prop.begin(), v_prop.end(), 0.);
      fill(n_prop.begin(), n_prop.end(), 0);
      for (auto &a: rgb_prop)
        fill(a.begin(), a.end(), 0.);

      // mapping correction
      int offset3 = static_cast<int>( round(0.5 * (nim->dim[3] - nim_tmp->dim[3])) );
      int offset2 = static_cast<int>( round(0.5 * (nim->dim[2] - nim_tmp->dim[2])) );
      int offset1 = static_cast<int>( round(0.5 * (nim->dim[1] - nim_tmp->dim[1])) );

      int c = 0;
      for (int k=0; k<nim_tmp->dim[3]; k++) {
        int kk = static_cast<int>( round(nim_tmp->pixdim[3] * k * conver_fac * vlen_1) );
        kk += offset3;

        for (int j=0; j<nim_tmp->dim[2]; j++) {
          int jj = static_cast<int>( round(nim_tmp->pixdim[2] * j * conver_fac * vlen_1) );
          jj += offset2;

          for (int i=0; i<nim_tmp->dim[1]; i++) {
            int ii = static_cast<int>( round(nim_tmp->pixdim[1] * i * conver_fac * vlen_1) );
            ii += offset1;

            tagint itag = find_tag(brn,ii,jj,kk);
            if (itag == -1) {
              //printf("Warning: a tag cannot be assigned for the voxels of the mri file. \n");
              c++;
              continue;
            }

            int vid = map[itag];

            double int16coef = 1.0 / (static_cast<double>(INT16_MAX)
                                      - static_cast<double>(INT16_MIN));
            double int32coef = 1.0 / (static_cast<double>(INT32_MAX)
                                      - static_cast<double>(INT32_MIN));
            // if it is in the partition
            if (vid != -1) {
              if (nim_tmp->datatype == DT_UINT8)
                v_prop[vid] += static_cast<double>(ptr8[c]);
              else if (nim_tmp->datatype == DT_INT16)
                v_prop[vid] += ( static_cast<double>(ptr16[c] - INT16_MIN) ) * int16coef;
              else if (nim_tmp->datatype == DT_INT32)
                v_prop[vid] += ( static_cast<double>(ptr32[c] - INT32_MIN) ) * int32coef;
              else if (nim_tmp->datatype == DT_FLOAT32)
                v_prop[vid] += static_cast<double>(ptrf[c]);
              else if (nim_tmp->datatype == DT_FLOAT64)
                v_prop[vid] += static_cast<double>(ptrd[c]);
              else if (nim_tmp->datatype == DT_RGB24) {
                rgb_prop[0][vid] += static_cast<double>(ptr_rgb[3*c]);
                rgb_prop[1][vid] += static_cast<double>(ptr_rgb[3*c + 1]);
                rgb_prop[2][vid] += static_cast<double>(ptr_rgb[3*c + 2]);
              }

              n_prop[vid]++;

              //printf("HERE: proc %i: vid = %i, v_prop[vid] = %g, ptr16[c] = %i, c = %i, %i \n",
                     //brn->me, vid, v_prop[vid], ptr16[c], c, nim_tmp->datatype);

            }

            c++;
          }
        }
      }

      // set voxel properties based on the nifti_image data
      for (int i=0; i<nall; i++) {

        if (n_prop[i] > 0) {
          double dum = 1.0 / n_prop[i];
          v_prop[i] *= dum;
          rgb_prop[0][i] *= dum;
          rgb_prop[1][i] *= dum;
          rgb_prop[2][i] *= dum;
        }

        double coef = 1.0 / (1.0 - thres_val);
        double coef_rgb = 1.0 / (max_val - thres_val);
        /* ----------------------------------------------------------------------
         * criteria based on mri file
         * setup all types and groups from a single file
         * ----------------------------------------------------------------------*/
        if (!mri_arg[tis][0].compare("all")) {
          if (v_prop[i] <= 0) {
            type[i] = tissue[EMP];
            for (int ag_id=0; ag_id<num_agents; ag_id++)
              agent[ag_id][i] = 0.0;
          }

          else if (v_prop[i] < thres_val) {
            type[i] = tissue[CSF];
            for (int ag_id=0; ag_id<num_agents; ag_id++)
              agent[ag_id][i] = 0.0;
          }

          else if (v_prop[i] > thres_val) {
            type[i] = tissue[GM];
            agent[neu][i] = (v_prop[i] - thres_val) * coef * max_val;
          }
        }
        /* ----------------------------------------------------------------------
         * setup WM type from wm file
         * ----------------------------------------------------------------------*/
        else if (!mri_arg[tis][0].compare("wm")) {
          if (v_prop[i] > thres_val) {
            type[i] = tissue[WM];
            agent[neu][i] += (v_prop[i] - thres_val) * coef * max_val;
          }
        }
        /* ----------------------------------------------------------------------
         * setup GM type from gm file
         * ----------------------------------------------------------------------*/
        else if (!mri_arg[tis][0].compare("gm")) {
          if (v_prop[i] > thres_val) {
            type[i] = tissue[GM];
            agent[neu][i] += (v_prop[i] - thres_val) * coef * max_val;
          }
        }
        /* ----------------------------------------------------------------------
         * setup CSF type from csf file
         * ----------------------------------------------------------------------*/
        else if (!mri_arg[tis][0].compare("csf")) {
          if (v_prop[i] > thres_val) {
            type[i] = tissue[CSF];
            for (int ag_id=0; ag_id<num_agents; ag_id++)
              agent[ag_id][i] = 0.0;
          }
        }
        /* ----------------------------------------------------------------------
         * setup groups from group file
         * ----------------------------------------------------------------------*/
        else if (!mri_arg[tis][0].compare("group")) {
          double fractpart, intpart;
          fractpart = modf (v_prop[i], &intpart);
          //printf("HERE0: %i %g \n",i,v_prop[i]);
          if (fractpart != 0) continue;
          group[i] = static_cast<int>( v_prop[i] );
        }
        /* ----------------------------------------------------------------------
         * setup diffusion tensor from RGB file
         * ----------------------------------------------------------------------*/
        else if (!mri_arg[tis][0].compare("rgb")) {
          for (int j=0; j<3; j++)
            if (rgb_prop[j][i] > thres_val)
              Dtau[j][i] = (rgb_prop[j][i] - thres_val) * coef_rgb * brn->Dtau_max;
        }

      }

      nifti_image_free(nim_tmp);

    }

  }

  // set all values to zero for EMT_type voxels
  for (int i=0; i<nall; i++) {
    if (type[i] & tissue[EMP]) {
      for (int ag_id=0; ag_id<num_agents; ag_id++)
        agent[ag_id][i] = 0.0;
      Dtau[0][i] = Dtau[1][i] = Dtau[2][i] = 0.0;
    }

    else if (type[i] & tissue[CSF])
      Dtau[0][i] = Dtau[1][i] = Dtau[2][i] = brn->Dtau_max;

    //printf("HERE %li max = %g , %g %g %g \n",
      //     i, brn->Dtau_max, Dtau[0][i],Dtau[1][i],Dtau[2][i]);

  }

}
