#include "scenario_connectome.h"

using namespace std;

/* ----------------------------------------------------------------------*/
ScenarioConnectome::ScenarioConnectome(int narg, char** arg, int rk, int np) {
  me = rk;
  nproc = np;

  MPI_Comm_split(MPI_COMM_WORLD,0,me,&world);

  reset();

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
void ScenarioConnectome::reset() {
  nvoxel = 0;
  nlocal = nghost = nall = 0;
  step = Nrun = 0;
  Nlog = 1000;

  dt = 0.0;
  nevery = -1;
  vlen = vlen_1 = vlen_2 = 0.0;
  vvol = vvol_1 = 0.0;

  init_val.clear();
  init_val.resize(num_agents,-1);

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

  nim = NULL;

  newton_flux = 1;

  // set tissue
  tissue.clear();
  tissue.resize(num_types);
  for (int i=0; i<num_types; i++)
    tissue[i] = 1 << i;

  input = new Input();
  init = new Init();
  comm = new Comm(this);
  output = new Output();
  region = new Region();

}

/* ----------------------------------------------------------------------*/
void ScenarioConnectome::allocations() {
  for (auto &a: x) {
    a.clear();
    a.resize(nall);
  }

  tag.clear();
  tag.resize(nall);

  type.clear();
  type.resize(nall);

  group.clear();
  group.resize(nall);

  is_loc.clear();
  is_loc.resize(nall);

  for (auto &a: agent) {
    a.clear();
    a.resize(nall);
  }

  for (auto &a: deriv) {
    a.clear();
    a.resize(nall);
  }

  for (auto &a: arr_prop.Dtau) {
    a.clear();
    a.resize(nall);
  }

  for (int ag_id=0; ag_id<num_agents; ag_id++) {
    if (init_val[ag_id] >= 0.0)
      for (int i=0; i<nall; i++)
        set_agent(ag_id,i,init_val[ag_id],0);
  }

  // set all voxel types as EMP and groups as 0
  fill(type.begin(),type.end(),tissue[EMP]);
  fill(group.begin(),group.end(),0);

  for (auto &a: arr_prop.Dtau)
    fill(a.begin(),a.end(),prop.Dtau_max);

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
          double dum = 0.5 * (arr_prop.Dtau[d][i] + arr_prop.Dtau[d][j]) * del_phr;
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

/* ----------------------------------------------------------------------*/
int ScenarioConnectome::set_property(string key, string val) {

  if (!key.compare("diff_sAb")) prop.diff_sAb = stof(val);
  else if (!key.compare("kp")) prop.kp = stof(val);
  else if (!key.compare("kn")) prop.kn = stof(val);
  else if (!key.compare("ds")) prop.ds = stof(val);
  else if (!key.compare("df")) prop.df = stof(val);
  else if (!key.compare("es")) prop.es = stof(val);
  else if (!key.compare("diff_mic")) prop.diff_mic = stof(val);
  else if (!key.compare("sens_s")) prop.sens_s = stof(val);
  else if (!key.compare("sens_f")) prop.sens_f = stof(val);
  else if (!key.compare("Ha")) prop.Ha = stof(val);
  else if (!key.compare("ka")) prop.ka = stof(val);
  else if (!key.compare("dnt")) prop.dnt = stof(val);
  else if (!key.compare("ktau")) prop.ktau = stof(val);
  else if (!key.compare("kphi")) prop.kphi = stof(val);
  else if (!key.compare("ephi")) prop.ephi = stof(val);
  else if (!key.compare("C_cir")) {
    prop.C_cir = stof(val);
    init_val[cir] = prop.C_cir;
  }
  else if (!key.compare("c_cir")) prop.c_cir = stof(val);
  else if (!key.compare("tau_cir")) prop.tau_cir = stof(val);
  else if (!key.compare("diff_tau")) prop.diff_tau = stof(val);
  else if (find_agent(key) >= 0) init_val[find_agent(key)] = stof(val);
  else return 0;

  return 1;
}

/* ----------------------------------------------------------------------*/
int ScenarioConnectome::find_agent(string str) {
  int ag_found = -1;

  for (int ag_id=0; ag_id<num_agents; ag_id++)
    if (!str.compare(ag_str[ag_id]))
      ag_found = ag_id;

  return ag_found;
}

/* ----------------------------------------------------------------------*/
void ScenarioConnectome::set_parameters() {
  prop.D_sAb = prop.diff_sAb * vlen_2;
  prop.D_mic = prop.diff_mic * vlen_2;
  prop.cs = prop.sens_s * vlen_2;
  prop.cf = prop.sens_f * vlen_2;
  prop.omega_cir = 2.0 * PI / prop.tau_cir;
  prop.Dtau_max = prop.diff_tau * vlen_2;
}

/* ----------------------------------------------------------------------
 * Define the system topology based on the mri nifti image (.nii)
 * ----------------------------------------------------------------------*/
void ScenarioConnectome::mri_topology(nifti_image *nim) {
  if (!nim)
    return;

  double conver_fac = 1.0;
  if (nim->xyz_units == NIFTI_UNITS_METER)
    conver_fac = 1.e6;
  else if (nim->xyz_units == NIFTI_UNITS_MM)
    conver_fac = 1.e3;
  else if (nim->xyz_units == NIFTI_UNITS_MICRON)
    conver_fac = 1.0;

  /* -------------------------------------------------------
   * set from restart
   * ------------------------------------------------------- */
  if (!init->mri_arg[0][0].compare("restart")) {
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

            tagint itag = init->find_tag(this,ii,jj,kk);
            if (itag == -1) {
              //printf("Warning: a tag cannot be assigned for the voxels of the mri file. \n");
              c++;
              continue;
            }

            int vid = init->map(itag);

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
        else if (find_agent(arg[i]) >= 0)
          set_agent(find_agent(arg[i]),j,v_prop[i][j],0);
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
    for (int tis=0; tis<init->mri_arg.size(); tis++) {
      nifti_image *nim_tmp = NULL;

      nim_tmp = nifti_image_read(init->mri_arg[tis][1].c_str(),1);
      thres_val = stof(init->mri_arg[tis][2]);
      max_val = stof(init->mri_arg[tis][3]);

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

            tagint itag = init->find_tag(this,ii,jj,kk);
            if (itag == -1) {
              //printf("Warning: a tag cannot be assigned for the voxels of the mri file. \n");
              c++;
              continue;
            }

            int vid = init->map(itag);

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
        if (!init->mri_arg[tis][0].compare("all")) {
          if (v_prop[i] <= 0) {
            type[i] = tissue[EMP];
            for (int ag_id=0; ag_id<num_agents; ag_id++)
              set_agent(ag_id,i,0.0,0);
          }

          else if (v_prop[i] < thres_val) {
            type[i] = tissue[CSF];
            for (int ag_id=0; ag_id<num_agents; ag_id++)
              set_agent(ag_id,i,0.0,0);
          }

          else if (v_prop[i] > thres_val) {
            type[i] = tissue[GM];
            int ag_id = find_agent("neu");
            set_agent(ag_id,i,(v_prop[i] - thres_val) * coef * max_val,0);
          }
        }
        /* ----------------------------------------------------------------------
         * setup WM type from wm file
         * ----------------------------------------------------------------------*/
        else if (!init->mri_arg[tis][0].compare("wm")) {
          if (v_prop[i] > thres_val) {
            type[i] = tissue[WM];
            int ag_id = find_agent("neu");
            set_agent(ag_id,i,(v_prop[i] - thres_val) * coef * max_val,1);
          }
        }
        /* ----------------------------------------------------------------------
         * setup GM type from gm file
         * ----------------------------------------------------------------------*/
        else if (!init->mri_arg[tis][0].compare("gm")) {
          if (v_prop[i] > thres_val) {
            type[i] = tissue[GM];
            int ag_id = find_agent("neu");
            set_agent(ag_id,i,(v_prop[i] - thres_val) * coef * max_val,1);
          }
        }
        /* ----------------------------------------------------------------------
         * setup CSF type from csf file
         * ----------------------------------------------------------------------*/
        else if (!init->mri_arg[tis][0].compare("csf")) {
          if (v_prop[i] > thres_val) {
            type[i] = tissue[CSF];
            for (int ag_id=0; ag_id<num_agents; ag_id++)
              set_agent(ag_id,i,0.0,0);
          }
        }
        /* ----------------------------------------------------------------------
         * setup groups from group file
         * ----------------------------------------------------------------------*/
        else if (!init->mri_arg[tis][0].compare("group")) {
          double fractpart, intpart;
          fractpart = modf (v_prop[i], &intpart);
          //printf("HERE0: %i %g \n",i,v_prop[i]);
          if (fractpart != 0) continue;
          group[i] = static_cast<int>( v_prop[i] );
        }
        /* ----------------------------------------------------------------------
         * setup diffusion tensor from RGB file
         * ----------------------------------------------------------------------*/
        else if (!init->mri_arg[tis][0].compare("rgb")) {
          for (int j=0; j<3; j++)
            if (rgb_prop[j][i] > thres_val)
              arr_prop.Dtau[j][i] = (rgb_prop[j][i] - thres_val) * coef_rgb * prop.Dtau_max;
        }

      }

      nifti_image_free(nim_tmp);

    }

  }

  // set all values to zero for EMP_type voxels
  for (int i=0; i<nall; i++) {
    if (type[i] & tissue[EMP]) {
      for (int ag_id=0; ag_id<num_agents; ag_id++)
        set_agent(ag_id,i,0.0,0);
      arr_prop.Dtau[0][i] = arr_prop.Dtau[1][i] = arr_prop.Dtau[2][i] = 0.0;
    }

    else if (type[i] & tissue[CSF])
      arr_prop.Dtau[0][i] = arr_prop.Dtau[1][i] = arr_prop.Dtau[2][i] = prop.Dtau_max;

  }

}

/* ----------------------------------------------------------------------*/
int ScenarioConnectome::dump_specific(vector<string> arg) {

  int dsize = 3;
  tagint c = 3;
  while (c < arg.size()) {
    if (!arg[c].compare("type"))
      dsize++;
    else if (!arg[c].compare("group"))
      dsize++;
    else if (!arg[c].compare("rgb")) {
      if (c != arg.size() - 1) {
        printf("Error: dump_mri: rgb keyword should be the last. \n");
        exit(1);
      }
      dsize += 3;
    }
    else if (!arg[c].compare("me"))
      dsize++;
    else if (find_agent(arg[c]) >= 0)
      dsize++;
    else
      return 0;
    c++;
  }

  vector<int> rcounts(nproc), displs(nproc);
  vector<double> send_buf(nlocal*dsize), recv_buf(nvoxel*dsize);

  // pack
  c = 0;
  for (int kk=1; kk<nvl[2]+1; kk++)
    for (int jj=1; jj<nvl[1]+1; jj++)
      for (int ii=1; ii<nvl[0]+1; ii++) {
        int i = find_id(ii,jj,kk);
        send_buf[c++] = x[0][i]; // x
        send_buf[c++] = x[1][i]; // y
        send_buf[c++] = x[2][i]; // z

        int aid = 3;
        while (aid < arg.size()) {
          int ag_id = find_agent(arg[aid]);
          if (!arg[aid].compare("type")) {
            int tis;
            for (tis=0; tis<num_types; tis++)
              if (type[i] & tissue[tis])
                break;
            send_buf[c++] = ubuf(tis).d;
          }
          else if (!arg[aid].compare("group"))
            send_buf[c++] = ubuf(group[i]).d;
          else if (!arg[aid].compare("rgb")) {
            send_buf[c++] = arr_prop.Dtau[0][i];
            send_buf[c++] = arr_prop.Dtau[1][i];
            send_buf[c++] = arr_prop.Dtau[2][i];
          }
          else if (!arg[aid].compare("me"))
            send_buf[c++] = ubuf(me).d;
          else if (ag_id >= 0)
            send_buf[c++] = get_agent(ag_id,i);
          aid++;
        }
      }

  MPI_Gather(&nlocal,1,MPI_INT,&rcounts[0],1,MPI_INT,0,world);

  if (!me) {
    int offset = 0;
    for (int i = 0; i < nproc; i++) {
      rcounts[i] *= dsize;
      displs[i] = offset;
      offset += rcounts[i];
    }
  }

  MPI_Gatherv(&send_buf[0],nlocal*dsize,MPI_DOUBLE,
              &recv_buf[0],&rcounts[0],&displs[0],MPI_DOUBLE,0,world);

  // unpack and print on the root
  if (!me) {
    const int dims5[] = {5, nv[0], nv[1], nv[2], 1, dsize-3, 1, 1};
    const int dims3[] = {3, nv[0], nv[1], nv[2], 1, 1, 1, 1};

    nifti_image *nim;

    if (dsize > 4)
      nim = output->nifti_image_setup(this,arg, dims5, NIFTI_INTENT_VECTOR);
    else
      nim = output->nifti_image_setup(this,arg, dims3, NIFTI_INTENT_NONE);

    float* data = (float*) nim->data;

    for (tagint i=0; i<nvoxel; i++) {
      tagint c = i * dsize;

      int ii = static_cast<int>( round((recv_buf[c++] - boxlo[0]) * vlen_1 - 0.5) );
      int jj = static_cast<int>( round((recv_buf[c++] - boxlo[1]) * vlen_1 - 0.5) );
      int kk = static_cast<int>( round((recv_buf[c++] - boxlo[2]) * vlen_1 - 0.5) );

      int aid = 3;
      while (aid < arg.size()) {

        tagint cnim = ii + nim->nx * ( jj + nim->ny * (kk + nim->nz * (aid-3) ) );
        int ag_id = find_agent(arg[aid]);

        if (!arg[aid].compare("type"))
          data[cnim] = (float) ubuf(recv_buf[c++]).i;
        else if (!arg[aid].compare("group"))
          data[cnim] = (float) ubuf(recv_buf[c++]).i;
        else if (!arg[aid].compare("rgb")) {
          data[cnim] = (float) recv_buf[c++];
          cnim = ii + nim->nx * ( jj + nim->ny * (kk + nim->nz * (aid-2) ) );
          data[cnim] = (float) recv_buf[c++];
          cnim = ii + nim->nx * ( jj + nim->ny * (kk + nim->nz * (aid-1) ) );
          data[cnim] = (float) recv_buf[c++];
        }
        else if (!arg[aid].compare("me"))
          data[cnim] = (float) ubuf(recv_buf[c++]).i;
        else if (ag_id >= 0)
          data[cnim] = (float) recv_buf[c++];

        aid++;
      }
    }

    nifti_image_write(nim);
    nifti_image_free(nim);

  }

  // successful mri output
  return 1;
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
