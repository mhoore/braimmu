#include "virtualbrain.h"
#include "input.h"

using namespace std;

/* ---------------------------------------------------------------------- */
Input::Input() {
}

/* ----------------------------------------------------------------------*/
Input::~Input() {
}

/* ----------------------------------------------------------------------*/
void Input::file(const char *fname, class VirtualBrain *brn) {
  string str;
  string comment ("#");

  fr.open(fname);

  do {
    getline(fr, str);

    if (!brn->me) {
      ofstream logfile;
      logfile.open (flog, ios::app);
      logfile << str << endl;
      logfile.close();
    }

    // ignore comments
    size_t found = str.find(comment);
    if (found != string::npos)
      str.erase(found,string::npos - found);
    if (str.size() == 0) continue;

    parse(str);

    execute_command(brn);

  } while (!fr.eof());

  fr.close();
}

/* ----------------------------------------------------------------------*/
void Input::parse(string str) {
  istringstream buf(str);

  buf >> command;

  arg.clear();
  for(string word; buf >> word;)
    arg.push_back(word);

  narg = arg.size();

}

/* ----------------------------------------------------------------------*/
void Input::execute_command(VirtualBrain *brn) {

  /* Command: timestep - the timestep of the simulation.
   * syntax: timestep dt
   * example: timestep 0.001
   * arguments: dt = timestep in (day) */
  if (!command.compare("timestep")) {
    if (narg < 1){
      printf("Error: timestep \n");
      exit(1);
    }

    brn->dt = stof(arg[0]);
  }

  /* Command: nevery - the period of long-term communications
   * for neural activity.
   * syntax: nevery nevery_value
   * example: nevery 100
   * arguments: nevery_value = integer multiples of dt */
  else if (!command.compare("nevery")) {
    if (narg < 1){
      printf("Error: nevery \n");
      exit(1);
    }

    brn->nevery = stoi(arg[0]);
  }

  /* Command: box - setting the simulation box boundaries.
   * syntax: box xlo xhi ylo yhi zlo zhi
   * example: box -1.0e4 1.0e4  -1.0e4 1.0e4  -1.0e4 1.0e4
   * arguments: low and high boundaries in (um) */
  else if (!command.compare("box")) {
    if (narg < 6){
      printf("Error: box \n");
      exit(1);
    }

    brn->boxlo[0] = stof(arg[0]);
    brn->boxhi[0] = stof(arg[1]);
    brn->boxlo[1] = stof(arg[2]);
    brn->boxhi[1] = stof(arg[3]);
    brn->boxlo[2] = stof(arg[4]);
    brn->boxhi[2] = stof(arg[5]);
  }

  /* Command: voxel_size - setting the lattice unit length.
   * syntax: voxel_size dx
   * example: voxel_size 1.0e3
   * arguments: voxel size in (um) */
  else if (!command.compare("voxel_size")) {
    if (narg < 1){
      printf("Error: voxel_size \n");
      exit(1);
    }

    brn->vlen = stof(arg[0]);
  }

  /* Command: partition - setting the partitioning for parallel
   * simulation, in each dimension, nx,ny,nz. The total number of
   * processors is the multiplication of them (nproc = nx*ny*nz).
   * syntax: partition nx ny nz
   * example: partition 4 2 2
   * arguments: minimum number is 1 */
  else if (!command.compare("partition")) {
    if (narg < 3){
      printf("Error: partition \n");
      exit(1);
    }

    brn->npart[0] = stoi(arg[0]);
    brn->npart[1] = stoi(arg[1]);
    brn->npart[2] = stoi(arg[2]);
  }

  else if (!command.compare("num_conn_max")) {
    if (narg < 1){
      printf("Error: num_conn_max \n");
      exit(1);
    }

    brn->num_conn_max = stoi(arg[0]);
  }

  /* Command: parameters - setting the initial values of the parameters.
   * syntax: parameters keyword value
   * example: parameters neu 6.e7 es 2.0
   * arguments: keywords can be any constant, or initial value of
   * agents in the model */
  else if (!command.compare("parameters")) {
    if (narg < 1) {
      printf("Error: parameters \n");
      exit(1);
    }

    if (!read_parameters(brn)) {
        printf("Error: read parameters \n");
        exit(1);
      }
  }

  /* Command: newton_flux: for calculating fluxes only once
   * between each pair of voxels.
   * syntax: newton_flux keyword
   * example: newton_flux no
   * arguments:
   * keywords: yes, no
   * default: yes */
  else if (!command.compare("newton_flux")) {
    if (narg < 1) {
      printf("Error: newton_flux \n");
      exit(1);
    }

    if (!arg[0].compare("yes")) brn->newton_flux = 1;
    else if (!arg[0].compare("no")) brn->newton_flux = 0;
    else {
      printf("Error: newton_flux keyword not recognized! \n");
      exit(1);
    }
  }

  /* Command: read_mri - assigning an mri nifti image file (.nii) for
   * the topology of the system.
   * syntax: read_mri keyword file.nii thres_val max_val
   * example: read_mri gm gm.nii 0.1 1.0e9
   * arguments:
   * keywords: gm, wm, csf, group, restart, all
   * thres_val: threshold value between 0 and 1
   * max_val: maximum number of neurons */
  else if (!command.compare("read_mri")) {
    if (narg < 4){
      printf("Error: read_mri \n");
      exit(1);
    }

    read_mri(brn);

  }

  /* Command: setup - initialization point the system.
   * syntax: setup
   * example:
   * arguments: setup must come after many commands.
   * After setup command, the system cannot be set any longer. */
  else if (!command.compare("setup")) {
    brn->init->setup(brn);
  }

  /* Command: region - setting a region
   * syntax: region region_style arguments in/out keyword value
   * example: region sphere 0.0 0.0 0.0 5.0e4 out type 1
   * arguments: region styles: sphere, block, etc.
   * simulation constants or agents initial values to be
   * set in the region. in/out means inside or outside the region. */
  else if (!command.compare("region")) {
    if (narg < 1){
      printf("Error: region \n");
      exit(1);
    }

    read_region(brn);

  }

  else if (!command.compare("restart")) {
    if (narg < 2) {
      printf("Error: restart \n");
      exit(1);
    }

    if (!read_restart(brn)) {
        printf("Error: read restart \n");
        exit(1);
      }
  }

  else if (!command.compare("statistics")) {
    if (narg < 2) {
      printf("Error: statistics \n");
      exit(1);
    }

    if (!read_statistics(brn)) {
        printf("Error: read statistics \n");
        exit(1);
      }
  }

  /* Command: dump - setting dump styles for outputing the results.
   * syntax: dump dump_style nevery filename keywords
   * example: dump mri 1000 dump.nii type neu
   * arguments: dump styles are "mri" and "txt" */
  else if (!command.compare("dump")) {
    if (narg < 3) {
      printf("Error: dump \n");
      exit(1);
    }

    if (!read_dump(brn)) {
        printf("Error: read dump \n");
        exit(1);
      }
  }

  else if (!command.compare("log")) {
    if (narg < 1){
      printf("Error: log \n");
      exit(1);
    }

    brn->Nlog = stoi(arg[0]);
  }

  else if (!command.compare("run")) {
    if (narg < 1){
      printf("Error: run \n");
      exit(1);
    }

    brn->Nrun = stoi(arg[0]);
  }

  else if (!command.compare("balance")) {
    if (narg < 2){
      printf("Error: balance \n");
      exit(1);
    }

    brn->comm->b_dim = arg[0];
    brn->comm->b_itr = stoi(arg[1]);
  }

  //else if (!command.compare("region"))
    //region();

  else {
    //printf("Error: Unknown command %s \n",command.c_str());
    printf("Error: Unknown command \n");
    exit(1);
  }

}

/* ----------------------------------------------------------------------*/
int Input::read_parameters(VirtualBrain *brn) {
  int c = 0;


  do {
    if (c+1 >= narg)
      return 0;
    if ( !brn->set_property(arg[c],arg[c+1]) )
      return 0;
    c += 2;
  } while (c < narg);

  return 1;
}

/* ----------------------------------------------------------------------
 * Read MRI image (.nii) to define the spatial domain of the simulation.
 * ----------------------------------------------------------------------*/
void Input::read_mri(VirtualBrain *brn) {
  brn->init->mri_arg.push_back(vector<string>());
  int i = brn->init->mri_arg.size() - 1;

  for (int j=0; j<narg; j++)
    brn->init->mri_arg[i].push_back(arg[j]);

}

/* ----------------------------------------------------------------------*/
void Input::read_region(VirtualBrain *brn) {
  brn->region->reg_arg.push_back(vector<string>());
  int i = brn->region->reg_arg.size() - 1;

  for (int j=0; j<narg; j++)
    brn->region->reg_arg[i].push_back(arg[j]);

}

/* ----------------------------------------------------------------------*/
int Input::read_restart(VirtualBrain *brn) {
  brn->output->revery = stoi(arg[0]);
  brn->output->rname = arg[1];

  // make a clean file
  if (!brn->me) {
    fstream fb;
    fb.open(brn->output->rname,ios::binary|ios::out|ios::trunc);
    fb.close();
  }

  return 1;
}

/* ----------------------------------------------------------------------*/
int Input::read_statistics(VirtualBrain *brn) {
  brn->output->severy = stoi(arg[0]);
  brn->output->sname = arg[1];

  // make a clean file
  if (!brn->me) {
    FILE *fw;
    fw = fopen(brn->output->sname.c_str(),"w");
    fprintf(fw,"# statistics over all voxels\n");
    fprintf(fw,"# step region ");
    for (int ag_id=0; ag_id<brn->get_num_agents(); ag_id++)
      fprintf(fw,"%s ", brn->get_ag_str(ag_id).c_str());
    fprintf(fw,"\n");
    fclose(fw);
  }

  return 1;
}

/* ----------------------------------------------------------------------*/
int Input::read_dump(VirtualBrain *brn) {
  brn->output->dump_arg.push_back(vector<string>());
  int i = brn->output->dump_arg.size() - 1;

  for (int j=0; j<narg; j++)
    brn->output->dump_arg[i].push_back(arg[j]);

  brn->output->do_dump = true;

  // make a clean file
  if ( !brn->me && !(brn->output->dump_arg[i][0].compare("txt")) ) {
    FILE *fw;
    fw = fopen(brn->output->dump_arg[i][2].c_str(),"w");
    fclose(fw);
  }

  return 1;
}
