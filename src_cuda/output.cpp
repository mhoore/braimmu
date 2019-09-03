#include "virtualbrain.h"
#include "output.h"

#include <climits>
#include <unistd.h> // lseek,read,write
#include <fcntl.h> // open etc.

using namespace std;

/* ---------------------------------------------------------------------- */
Output::Output() {
  do_dump = false;
  severy = -1;
}

/* ----------------------------------------------------------------------*/
Output::~Output() {
  for (int i=0; i<dump_arg.size(); i++)
    dump_arg[i].clear();

  dump_arg.clear();
}

/* ----------------------------------------------------------------------*/
void Output::lammpstrj(VirtualBrain *brn) {

  MPI_Comm world = brn->world;

  int me = brn->me;
  int nproc = brn->nproc;

  int nlocal = brn->nlocal;
  tagint nvoxel = brn->nvoxel;

  auto &nvl = brn->nvl;

  auto &tag = brn->tag;
  auto &tissue = brn->tissue;
  auto &type = brn->type;
  auto &group = brn->group;

  auto &x = brn->x;

  int num_agents = brn->get_num_agents();

  int dsize = num_agents + 6; // NOTE: the number of data packed to the buffer
  vector<int> rcounts(nproc), displs(nproc);
  vector<double> send_buf(nlocal*dsize), recv_buf(nvoxel*dsize);

  // pack
  tagint c = 0;
  for (int kk=1; kk<nvl[2]+1; kk++)
    for (int jj=1; jj<nvl[1]+1; jj++)
      for (int ii=1; ii<nvl[0]+1; ii++) {
        int i = brn->find_id(ii,jj,kk);
        send_buf[c++] = ubuf(tag[i]).d;
        // type
        int tis;
        for (tis=0; tis<num_types; tis++)
          if (type[i] & tissue[tis])
            break;
        send_buf[c++] = ubuf(tis).d;
        send_buf[c++] = ubuf(group[i]).d;

        send_buf[c++] = x[0][i];
        send_buf[c++] = x[1][i];
        send_buf[c++] = x[2][i];

        for (int ag_id=0; ag_id<num_agents; ag_id++)
          send_buf[c++] = brn->get_agent(ag_id,i);
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

  // unpack and print
  if (!me) {
    fw = fopen("init.lammpstrj","w");

    fprintf(fw,"ITEM: TIMESTEP \n");
    fprintf(fw,"%i \n", brn->step);

    fprintf(fw,"ITEM: NUMBER OF ATOMS \n");
    fprintf(fw,TAGINT_FORMAT " \n", nvoxel);

    fprintf(fw,"ITEM: BOX BOUNDS pp pp pp \n");
    fprintf(fw,"%g %g \n", brn->boxlo[0],brn->boxhi[0]);
    fprintf(fw,"%g %g \n", brn->boxlo[1],brn->boxhi[1]);
    fprintf(fw,"%g %g \n", brn->boxlo[2],brn->boxhi[2]);

    fprintf(fw,"ITEM: ATOMS id type x y z ");
    for (int ag_id=0; ag_id<num_agents; ag_id++)
      fprintf(fw,"%s ", brn->get_ag_str(ag_id).c_str());
    fprintf(fw,"\n");

    tagint c = 0;
    for (int i=0; i<nvoxel; i++) {
      fprintf(fw,TAGINT_FORMAT " ", (tagint) ubuf(recv_buf[c++]).i); // tag
      fprintf(fw,"%i ", (int) ubuf(recv_buf[c++]).i); // type
      c++;
      //fprintf(fw,"%i ", (int) ubuf(recv_buf[c++]).i); // group
      fprintf(fw,"%g ", recv_buf[c++]); // x
      fprintf(fw,"%g ", recv_buf[c++]); // y
      fprintf(fw,"%g ", recv_buf[c++]); // z
      for (int ag_id=0; ag_id<num_agents; ag_id++)
        fprintf(fw,"%.3f ", recv_buf[c++]); // agents
      fprintf(fw,"\n");
    }

    fclose(fw);

  }

}

/* ----------------------------------------------------------------------*/
void Output::restart(VirtualBrain *brn) {
  if (brn->step % revery != 0) return;

  MPI_Comm world = brn->world;

  int me = brn->me;
  int nproc = brn->nproc;

  int nlocal = brn->nlocal;
  tagint nvoxel = brn->nvoxel;

  auto &tag = brn->tag;
  auto &tissue = brn->tissue;
  auto &type = brn->type;
  auto &group = brn->group;

  auto &nvl = brn->nvl;

  int num_agents = brn->get_num_agents();

  int dsize = num_agents + 3;

  vector<int> rcounts, displs;
  vector<double> send_buf, recv_buf;

  send_buf.clear();
  send_buf.resize(nlocal*dsize);

  // pack
  tagint c = 0;
  for (int kk=1; kk<nvl[2]+1; kk++)
    for (int jj=1; jj<nvl[1]+1; jj++)
      for (int ii=1; ii<nvl[0]+1; ii++) {
        int i = brn->find_id(ii,jj,kk);
        send_buf[c++] = ubuf(tag[i]).d;

        // type
        int tis;
        for (tis=0; tis<num_types; tis++)
          if (type[i] & tissue[tis])
            break;
        send_buf[c++] = ubuf(tis).d;

        send_buf[c++] = ubuf(group[i]).d;

        for (int ag_id=0; ag_id<num_agents; ag_id++)
          send_buf[c++] = brn->get_agent(ag_id,i);
      }

  if (!me) {
    rcounts.clear();
    displs.clear();
    recv_buf.clear();

    rcounts.resize(nproc);
    displs.resize(nproc);
    recv_buf.resize(nvoxel*dsize);
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

  // unpack and print
  if (!me) {
    tagint buft;
    int bufi;
    double bufd;

    //write
    fb.open(rname,ios::binary|ios::out|ios::trunc);
    fb.write((char *)&brn->step,sizeof(brn->step));
    fb.write((char *)&nvoxel,sizeof(nvoxel));

    fb.write((char *)&brn->vlen,sizeof(brn->vlen));

    fb.write((char *)&brn->boxlo[0],sizeof(brn->boxlo[0]));
    fb.write((char *)&brn->boxhi[0],sizeof(brn->boxhi[0]));

    fb.write((char *)&brn->boxlo[1],sizeof(brn->boxlo[1]));
    fb.write((char *)&brn->boxhi[1],sizeof(brn->boxhi[1]));

    fb.write((char *)&brn->boxlo[2],sizeof(brn->boxlo[2]));
    fb.write((char *)&brn->boxhi[2],sizeof(brn->boxhi[2]));

    tagint c = 0;
    for (int i=0; i<nvoxel; i++) {
      buft = (tagint) ubuf(recv_buf[c++]).i; // tag
      fb.write((char *)&buft,sizeof(buft));

      // type
      bufi = (int) ubuf(recv_buf[c++]).i;
      fb.write((char *)&bufi,sizeof(bufi)); // type

      bufi = (int) ubuf(recv_buf[c++]).i; // group
      fb.write((char *)&bufi,sizeof(bufi));

      for (int ag_id=0; ag_id<num_agents; ag_id++) {
        bufd = recv_buf[c++];
        fb.write((char *)&bufd,sizeof(bufd));
      }
    }

    fb.close();

    /*
    //read and check
    fw = fopen("test.txt","w");
    fb.open(dname,ios::binary|ios::in|ios::app);
    fb.seekg(0);
    fb.read((char *)&bufi,sizeof(bufi)); // step
    fprintf(fw,"%i \n", bufi);
    fb.read((char *)&buft,sizeof(buft)); //nvoxel
    fprintf(fw,TAGINT_FORMAT " \n", buft);
    for (i=0; i<nvoxel; i++) {
      fb.read((char *)&buft,sizeof(buft)); // tag
      fprintf(fw,TAGINT_FORMAT " ", buft);
      fb.read((char *)&bufi,sizeof(bufi)); // type
      fprintf(fw,"%i ", bufi);
      fb.read((char *)&bufd,sizeof(bufd)); // fAb
      fprintf(fw,"%g \n", bufd);
    }
    fb.close();
    fclose(fw);
    */

  }

}

/* ----------------------------------------------------------------------*/
void Output::statistics(VirtualBrain *brn) {
  if (brn->step % severy != 0) return;

  MPI_Comm world = brn->world;

  int me = brn->me;

  auto &tissue = brn->tissue;
  auto &type = brn->type;

  auto &nvl = brn->nvl;

  int num_agents = brn->get_num_agents();

  int nr = 2; // number of regions: parenchyma and CSF

  vector<vector<double>> agent_val(num_agents);
  for (auto &a: agent_val) {
    a.clear();
    a.resize(nr);
    fill(a.begin(), a.end(), 0.0);
  }

  vector<tagint> agent_num(nr,0);

  for (int kk=1; kk<nvl[2]+1; kk++)
    for (int jj=1; jj<nvl[1]+1; jj++)
      for (int ii=1; ii<nvl[0]+1; ii++) {
        int i = brn->find_id(ii,jj,kk);

        int c;
        if (type[i] & tissue[EMP]) continue;
        else if (type[i] & tissue[WM] || type[i] & tissue[GM])
          c = 0;
        else if (type[i] & tissue[CSF])
          c = 1;

        for (int ag_id=0; ag_id<num_agents; ag_id++)
          agent_val[ag_id][c] += brn->get_agent(ag_id,i);
        agent_num[c]++;
      }

  int size = nr*(num_agents+1);
  vector<double> send_buf(size);

  // pack
  tagint c = 0;
  for (int i=0; i<nr; i++) {
    send_buf[c++] = ubuf(agent_num[i]).d;
    for (int ag_id=0; ag_id<num_agents; ag_id++)
      send_buf[c++] = agent_val[ag_id][i];
  }

  vector<double> recv_buf(size);

  MPI_Reduce(&send_buf[0],&recv_buf[0],size,MPI_DOUBLE,MPI_SUM,0,world);

  if (!me) {
    // unpack
    tagint c = 0;
    for (int i=0; i<nr; i++) {
      agent_num[i] = (tagint) ubuf(recv_buf[c++]).i;
      for (int ag_id=0; ag_id<num_agents; ag_id++)
        agent_val[ag_id][i] = recv_buf[c++] / agent_num[i];
    }

    // print
    fw = fopen(sname.c_str(),"a");

    // PAR
    fprintf(fw,"%i ", brn->step);
    fprintf(fw,"PAR ");
    for (int ag_id=0; ag_id<num_agents; ag_id++)
      fprintf(fw,"%g ", agent_val[ag_id][0]); // agent_values
    fprintf(fw,"\n");

    // CSF
    fprintf(fw,"%i ", brn->step);
    fprintf(fw,"CSF ");
    for (int ag_id=0; ag_id<num_agents; ag_id++)
      fprintf(fw,"%g ", agent_val[ag_id][1]); // agent_values
    fprintf(fw,"\n");

    fclose(fw);

  }

}

/* ----------------------------------------------------------------------*/
void Output::statistics_sphere(VirtualBrain *brn) {
  if (brn->step % severy != 0) return;

  MPI_Comm world = brn->world;

  int me = brn->me;

  double vlen = brn->vlen;
  double vlen_1 = brn->vlen_1;

  auto &x = brn->x;

  auto &nvl = brn->nvl;

  tagint nr = (tagint)( brn->lbox[0] * vlen_1 / 2 );

  int num_agents = brn->get_num_agents();

  vector<vector<double>> agent_val(num_agents);
  for (auto &a: agent_val) {
    a.clear();
    a.resize(nr);
    fill(a.begin(), a.end(), 0.0);
  }

  vector<tagint> agent_num(nr,0);

  for (int kk=1; kk<nvl[2]+1; kk++)
    for (int jj=1; jj<nvl[1]+1; jj++)
      for (int ii=1; ii<nvl[0]+1; ii++) {
        int i = brn->find_id(ii,jj,kk);
        int c = (tagint) (vlen_1 * sqrt(x[0][i] * x[0][i]
                                      + x[1][i] * x[1][i]
                                      + x[2][i] * x[2][i]));
        //c = (tagint) (vlen_1 * (x[0][i] - brn->boxlo[0]));
        if (c >= nr) continue;

        for (int ag_id=0; ag_id<num_agents; ag_id++)
          agent_val[ag_id][c] += brn->get_agent(ag_id,i);
        agent_num[c]++;
      }

  int size = nr*(num_agents+1);
  vector<double> send_buf(size);

  // pack
  tagint c = 0;
  for (int i=0; i<nr; i++) {
    send_buf[c++] = ubuf(agent_num[i]).d;
    for (int ag_id=0; ag_id<num_agents; ag_id++)
      send_buf[c++] = agent_val[ag_id][i];
  }

  vector<double> recv_buf(size);

  MPI_Reduce(&send_buf[0],&recv_buf[0],size,MPI_DOUBLE,MPI_SUM,0,world);

  if (!me) {
    // unpack
    tagint c = 0;
    for (int i=0; i<nr; i++) {
      agent_num[i] = (tagint) ubuf(recv_buf[c++]).i;
      for (int ag_id=0; ag_id<num_agents; ag_id++)
        agent_val[ag_id][i] = recv_buf[c++] / agent_num[i];
    }

    // print
    fw = fopen(sname.c_str(),"a");

    fprintf(fw,"ITEM: TIMESTEP \n");
    fprintf(fw,"%i \n", brn->step);

    fprintf(fw,"ITEM: NUMBER OF BINS \n");
    fprintf(fw,TAGINT_FORMAT " \n", nr);

    fprintf(fw,"ITEM: BIN id r ");
    for (int ag_id=0; ag_id<num_agents; ag_id++)
      fprintf(fw,"%s ", brn->get_ag_str(ag_id).c_str());
    fprintf(fw,"\n");

    for (int i=0; i<nr; i++) {
      fprintf(fw,TAGINT_FORMAT " ", i); // id
      fprintf(fw,"%g ", vlen * i); // r
      for (int ag_id=0; ag_id<num_agents; ag_id++)
        fprintf(fw,"%g ", agent_val[ag_id][i]); // agent_values
      fprintf(fw,"\n");
    }

    fclose(fw);

  }

}

/* ----------------------------------------------------------------------*/
void Output::dump(VirtualBrain *brn) {
  for (int i=0; i<dump_arg.size(); i++) {

    if (brn->step % stoi(dump_arg[i][1]) != 0) continue;

    if (!dump_arg[i][0].compare("mri"))
      dump_mri(brn,dump_arg[i]);
    else if (!dump_arg[i][0].compare("txt"))
      dump_txt(brn,dump_arg[i]);
    else {
      printf("Error: dump %s unknown. \n", dump_arg[i][0].c_str());
      exit(1);
    }

  }

}

/* ----------------------------------------------------------------------*/
void Output::dump_txt(VirtualBrain *brn, vector<string> arg) {

  MPI_Comm world = brn->world;

  int me = brn->me;
  int nproc = brn->nproc;

  int nlocal = brn->nlocal;
  tagint nvoxel = brn->nvoxel;

  auto &x = brn->x;
  auto &tissue = brn->tissue;
  auto &type = brn->type;
  auto &group = brn->group;

  auto &nvl = brn->nvl;

  int dsize = 3;
  tagint c = 3;
  while (c < arg.size()) {
    if (!arg[c].compare("type"))
      dsize++;
    else if (!arg[c].compare("group"))
      dsize++;
    else if (brn->find_agent(arg[c]) >= 0)
      dsize++;
    else {
      printf("Error: dump argument \"%s\" unknown. \n", arg[c].c_str());
      exit(1);
    }
    c++;
  }

  vector<int> rcounts(nproc), displs(nproc);
  vector<double> send_buf(nlocal*dsize), recv_buf(nvoxel*dsize);

  // pack
  c = 0;
  for (int kk=1; kk<nvl[2]+1; kk++)
    for (int jj=1; jj<nvl[1]+1; jj++)
      for (int ii=1; ii<nvl[0]+1; ii++) {
        int i = brn->find_id(ii,jj,kk);
        send_buf[c++] = x[0][i]; // x
        send_buf[c++] = x[1][i]; // y
        send_buf[c++] = x[2][i]; // z

        int aid = 3;
        while (aid < arg.size()) {
          int ag_id = brn->find_agent(arg[aid]);
          if (!arg[aid].compare("type")) {
            // type
            int tis;
            for (tis=0; tis<num_types; tis++)
              if (type[i] & tissue[tis])
                break;
            send_buf[c++] = ubuf(tis).d;
          }
          else if (!arg[aid].compare("group"))
            send_buf[c++] = ubuf(group[i]).d;
          else if (ag_id >= 0)
            send_buf[c++] = brn->get_agent(ag_id,i);
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
    fw = fopen(arg[2].c_str(),"a");

    fprintf(fw,"ITEM: TIMESTEP \n");

    fprintf(fw,"%i \n", brn->step);

    fprintf(fw,"ITEM: NUMBER OF VOXELS \n");
    fprintf(fw,TAGINT_FORMAT " \n", nvoxel);

    fprintf(fw,"ITEM: VOXEL LENGTH \n");
    fprintf(fw,"%g \n", brn->vlen);

    fprintf(fw,"ITEM: BOUNDARIES \n");
    fprintf(fw,"%g %g \n", brn->boxlo[0], brn->boxhi[0]);
    fprintf(fw,"%g %g \n", brn->boxlo[1], brn->boxhi[1]);
    fprintf(fw,"%g %g \n", brn->boxlo[2], brn->boxhi[2]);

    fprintf(fw,"ITEM: VOXELS id ");

    int aid = 3;
    while (aid < arg.size()) {
      int ag_id = brn->find_agent(arg[aid]);
      if (!arg[aid].compare("type"))
        fprintf(fw,"type ");
      else if (!arg[aid].compare("group"))
        fprintf(fw,"group ");
      else if (ag_id >= 0)
        fprintf(fw,"%s ", brn->get_ag_str(ag_id).c_str());
      aid++;
    }
    fprintf(fw,"\n");

    for (int i=0; i<nvoxel; i++) {
      tagint c = i * dsize;
      tagint itag = find_tag(brn,recv_buf[c++],recv_buf[c++],recv_buf[c++]);

      fprintf(fw,TAGINT_FORMAT " ", itag); // tag

      int aid = 3;
      while (aid < arg.size()) {
        int ag_id = brn->find_agent(arg[aid]);
        if (!arg[aid].compare("type"))
          fprintf(fw,"%i ", (int) ubuf(recv_buf[c++]).i); // type
        else if (!arg[aid].compare("group"))
          fprintf(fw,"%i ", (int) ubuf(recv_buf[c++]).i); // group
        else if (ag_id >= 0)
          fprintf(fw,"%g ", recv_buf[c++]); // agent
        aid++;
      }
      fprintf(fw,"\n");

    }

    fclose(fw);

  }

}

/* ----------------------------------------------------------------------*/
void Output::dump_mri(VirtualBrain *brn, vector<string> arg) {

  /// check the arguments
  tagint c = 3;
  while (c < arg.size()) {
    if ( !arg[c].compare("type")
         || !arg[c].compare("group")
         || !arg[c].compare("me")
         || brn->find_agent(arg[c]) >= 0 )
      c++;
    else if (brn->dump_specific(arg))
        return;
    else {
      printf("Error: dump argument \"%s\" unknown. \n", arg[c].c_str());
      exit(1);
    }
  }

  MPI_Comm world = brn->world;

  int me = brn->me;
  int nproc = brn->nproc;

  int nlocal = brn->nlocal;
  tagint nvoxel = brn->nvoxel;

  auto &x = brn->x;
  auto &tissue = brn->tissue;
  auto &type = brn->type;
  auto &group = brn->group;

  auto &nvl = brn->nvl;

  int dsize = 3;
  c = 3;
  while (c < arg.size()) {
    if (!arg[c].compare("type"))
      dsize++;
    else if (!arg[c].compare("group"))
      dsize++;
    else if (!arg[c].compare("me"))
      dsize++;
    else if (brn->find_agent(arg[c]) >= 0)
      dsize++;
    c++;
  }

  vector<int> rcounts(nproc), displs(nproc);
  vector<double> send_buf(nlocal*dsize), recv_buf(nvoxel*dsize);

  // pack
  c = 0;
  for (int kk=1; kk<nvl[2]+1; kk++)
    for (int jj=1; jj<nvl[1]+1; jj++)
      for (int ii=1; ii<nvl[0]+1; ii++) {
        int i = brn->find_id(ii,jj,kk);
        send_buf[c++] = x[0][i]; // x
        send_buf[c++] = x[1][i]; // y
        send_buf[c++] = x[2][i]; // z

        int aid = 3;
        while (aid < arg.size()) {
          int ag_id = brn->find_agent(arg[aid]);
          if (!arg[aid].compare("type")) {
            int tis;
            for (tis=0; tis<num_types; tis++)
              if (type[i] & tissue[tis])
                break;
            send_buf[c++] = ubuf(tis).d;
          }
          else if (!arg[aid].compare("group"))
            send_buf[c++] = ubuf(group[i]).d;
          else if (!arg[aid].compare("me"))
            send_buf[c++] = ubuf(me).d;
          else if (ag_id >= 0)
            send_buf[c++] = brn->get_agent(ag_id,i);
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
    const int dims5[] = {5, brn->nv[0], brn->nv[1], brn->nv[2], 1, dsize-3, 1, 1};
    const int dims3[] = {3, brn->nv[0], brn->nv[1], brn->nv[2], 1, 1, 1, 1};

    nifti_image *nim;

    if (dsize > 4)
      nim = nifti_image_setup(brn,arg, dims5, NIFTI_INTENT_VECTOR);
    else
      nim = nifti_image_setup(brn,arg, dims3, NIFTI_INTENT_NONE);

    float* data = (float*) nim->data;

    auto &boxlo = brn->boxlo;
    double vlen_1 = brn->vlen_1;

    for (tagint i=0; i<brn->nvoxel; i++) {
      tagint c = i * dsize;

      int ii = static_cast<int>( round((recv_buf[c++] - boxlo[0]) * vlen_1 - 0.5) );
      int jj = static_cast<int>( round((recv_buf[c++] - boxlo[1]) * vlen_1 - 0.5) );
      int kk = static_cast<int>( round((recv_buf[c++] - boxlo[2]) * vlen_1 - 0.5) );

      int aid = 3;
      while (aid < arg.size()) {

        tagint cnim = ii + nim->nx * ( jj + nim->ny * (kk + nim->nz * (aid-3) ) );
        int ag_id = brn->find_agent(arg[aid]);

        if (!arg[aid].compare("type"))
          data[cnim] = (float) ubuf(recv_buf[c++]).i;
        else if (!arg[aid].compare("group"))
          data[cnim] = (float) ubuf(recv_buf[c++]).i;
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

}

/* ----------------------------------------------------------------------*/
void Output::sort_tag(VirtualBrain *brn, double *data, int dsize) {
  tagint i = 0;
  while (i < brn->nvoxel) {
    tagint ii = i*dsize;

    double xx[3];
    xx[0] = data[ii];
    xx[1] = data[ii+1];
    xx[2] = data[ii+2];

    tagint itag = find_tag(brn, xx[0], xx[1], xx[2]);

    if (i == itag)
      i++;
    else {
      // exchange data
      tagint it = itag*dsize;
      for (int j=0; j<dsize; j++) {
        double dum = data[it + j];
        data[it + j] = data[ii + j];
        data[ii + j] = dum;
      }
    }

  }

}

/* ----------------------------------------------------------------------
 * Find the tag of a voxel from its global location x, y, z
 * ----------------------------------------------------------------------*/
tagint Output::find_tag(VirtualBrain *brn, double x, double y, double z) {
  auto &nv = brn->nv;
  auto &boxlo = brn->boxlo;
  double vlen_1 = brn->vlen_1;

  int i = static_cast<int>( round((x - boxlo[0]) * vlen_1 - 0.5) );
  int j = static_cast<int>( round((y - boxlo[1]) * vlen_1 - 0.5) );
  int k = static_cast<int>( round((z - boxlo[2]) * vlen_1 - 0.5) );

  return i + nv[0] * (j + nv[1] * k);
}

/* ----------------------------------------------------------------------*/
nifti_image* Output::nifti_image_setup(VirtualBrain *brn, vector<string> arg,
                                       const int dims[], int intent) {
  nifti_image *nim;
  nim = nifti_make_new_nim(dims,DT_FLOAT32,1);

  // header and image file names (.nii)
  string str = "";
  str.append(arg[2]);
  str.append(to_string(brn->step));
  const char *prefix = str.c_str();

  if( nifti_set_filenames(nim, prefix, 1, 1) ) {
    printf("Error: nifti_set_filenames cannot set the names for %s . \n", prefix);
    exit(1);
  }

  for (int i=0; i<8; i++)
    nim->dim[i] = dims[i];

  nim->intent_p1 = nim->intent_p2 = nim->intent_p3 = 0.0;

  nim->intent_code = intent; // NIFTI_INTENT_VECTOR;

  nim->datatype = DT_FLOAT32;
  nim->nbyper = 4;

  nim->slice_start = 0;
  nim->pixdim[0] = 1;
  nim->pixdim[1] = brn->vlen;
  nim->pixdim[2] = brn->vlen;
  nim->pixdim[3] = brn->vlen;
  nim->pixdim[4] = brn->dt * stoi(arg[1]);
  nim->pixdim[5] = 1;
  nim->pixdim[6] = 1;
  nim->pixdim[7] = 1;

  nim->ndim = nim->dim[0];
  nim->nx = nim->dim[1];
  nim->ny = nim->dim[2];
  nim->nz = nim->dim[3];
  nim->nt = nim->dim[4];
  nim->nu = nim->dim[5];
  nim->nv = nim->dim[6];
  nim->nw = nim->dim[7];

  nim->dx = nim->pixdim[1];
  nim->dy = nim->pixdim[2];
  nim->dz = nim->pixdim[3];
  nim->dt = nim->pixdim[4];
  nim->du = nim->pixdim[5];
  nim->dv = nim->pixdim[6];
  nim->dw = nim->pixdim[7];

  nim->nvox = 1;
  for (int i=1; i<8; i++)
    nim->nvox *= nim->dim[i];

  nim->scl_slope = 0.0;
  nim->scl_inter = 0.0;

  nim->slice_end = 0.0;
  nim->slice_code = 0.0;

  nim->cal_min = 0.0;
  nim->cal_max = 0.0;

  nim->xyz_units = NIFTI_UNITS_MICRON;
  nim->time_units = NIFTI_UNITS_UNKNOWN;

  nim->slice_duration = 0.0;
  // initial time of the image
  nim->toffset = brn->step * brn->dt;

  str = "";
  for (int i=3; i<arg.size(); i++) {
    str.append(arg[i]);
    str.append(" ");
  }

  strcpy(nim->descrip, str.c_str());

  memset(nim->intent_name, 0, sizeof(nim->intent_name));
  memset(nim->aux_file, 0, sizeof(nim->aux_file));

  // voxel spatial position mapping scheme
  nim->qform_code = 0; // Method 2: using quaternions
  nim->sform_code = 0; // Method 3: using general affine transformation

  // quaternions - Method 2 of mapping
  nim->quatern_b = nim->quatern_c = nim->quatern_d = 0.0;
  nim->qoffset_x = nim->qoffset_y = nim->qoffset_z = 0.0;
  nim->qfac = 1.0;

  nim->cal_min = 0.0;
  nim->cal_max = 0.0;

  // not so important parameters
  nim->freq_dim = 0;
  nim->phase_dim = 0;
  nim->slice_dim = 0;
  nim->slice_code = 0;
  nim->slice_start = 0;
  nim->slice_end = 0;
  nim->slice_duration = 0.0;

  // method3 mapping
  nim->sto_xyz.m[0][0] = 1.0;
  nim->sto_xyz.m[0][1] = 0.0;
  nim->sto_xyz.m[0][2] = 0.0;
  nim->sto_xyz.m[0][3] = 0.0;
  nim->sto_xyz.m[1][0] = 0.0;
  nim->sto_xyz.m[1][1] = 1.0;
  nim->sto_xyz.m[1][2] = 0.0;
  nim->sto_xyz.m[1][3] = 0.0;
  nim->sto_xyz.m[2][0] = 0.0;
  nim->sto_xyz.m[2][1] = 0.0;
  nim->sto_xyz.m[2][2] = 1.0;
  nim->sto_xyz.m[2][3] = 0.0;

  nim->nifti_type = 1; // NIFTI-1 (1 file)

  return nim;

}
