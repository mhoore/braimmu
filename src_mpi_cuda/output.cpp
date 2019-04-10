#include <mpi.h>
#include "math.h"

#include "pointers.h"
#include "brain.h"
#include "output.h"

#include <climits>

#include <unistd.h> // lseek,read,write
#include <fcntl.h> // open etc.

using namespace std;
using namespace brain_NS;

/* ---------------------------------------------------------------------- */
Output::Output() {
  do_dump = false;
  severy = -1;
}

/* ----------------------------------------------------------------------*/
Output::~Output() {
  int i;

  for (i=0; i<dump_arg.size(); i++) {
    dump_arg[i].clear();
  }

  dump_arg.clear();
}

/* ----------------------------------------------------------------------*/
void Output::lammpstrj(Brain *brn) {
  tagint i,c;
  int dsize, ag_id;

  int *rcounts = NULL;
  int *displs = NULL;
  double *send_buf = NULL;
  double *recv_buf = NULL;

  MPI_Comm world = brn->world;

  int me = brn->me;
  int nproc = brn->nproc;

  int nlocal = brn->nlocal;
  int nall = brn->nall;
  tagint nvoxel = brn->nvoxel;

  tagint *tag = brn->tag;

  double **x = brn->x;

  dsize = num_agents + 5; // NOTE: the number of data packed to the buffer
  create(send_buf,nlocal*dsize,"send_buf");

  // pack
  c = 0;
  for (i=0; i<nall; i++) {
    if (!brn->is_loc[i]) continue; 
    send_buf[c++] = ubuf(tag[i]).d;
    send_buf[c++] = ubuf(brn->type[i]).d;

    send_buf[c++] = x[i][0];
    send_buf[c++] = x[i][1];
    send_buf[c++] = x[i][2];

    for (ag_id=0; ag_id<num_agents; ag_id++)
      send_buf[c++] = brn->agent[ag_id][i];
  }

  create(rcounts,nproc,"rcounts");
  create(displs,nproc,"displs");
  create(recv_buf,nvoxel*dsize,"recv_buf");

  MPI_Gather(&nlocal,1,MPI_INT,rcounts,1,MPI_INT,0,world);

  if (!me) {
    int offset = 0;
    for (i = 0; i < nproc; i++) {
      rcounts[i] *= dsize;
      displs[i] = offset;
      offset += rcounts[i];
    }
  }

  MPI_Gatherv(send_buf,nlocal*dsize,MPI_DOUBLE,recv_buf,rcounts,displs,MPI_DOUBLE,0,world);

  destroy(send_buf);

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
    for (ag_id=0; ag_id<num_agents; ag_id++)
      fprintf(fw,"%s ", ag_str[ag_id].c_str());
    fprintf(fw,"\n");

    c = 0;
    for (i=0; i<nvoxel; i++) {
      fprintf(fw,TAGINT_FORMAT " ", (tagint) ubuf(recv_buf[c++]).i); // tag
      fprintf(fw,"%i ", (int) ubuf(recv_buf[c++]).i); // type
      fprintf(fw,"%g ", recv_buf[c++]); // x
      fprintf(fw,"%g ", recv_buf[c++]); // y
      fprintf(fw,"%g ", recv_buf[c++]); // z
      for (ag_id=0; ag_id<num_agents; ag_id++)
        fprintf(fw,"%.3f ", recv_buf[c++]); // agents
      fprintf(fw,"\n");
    }

    fclose(fw);

  }

  destroy(rcounts);
  destroy(displs);
  destroy(recv_buf);

}

/* ----------------------------------------------------------------------*/
void Output::restart(Brain *brn) {
  if (brn->step % revery != 0) return;

  tagint i,c;
  int dsize, ag_id;

  int *rcounts = NULL;
  int *displs = NULL;
  double *send_buf = NULL;
  double *recv_buf = NULL;

  MPI_Comm world = brn->world;

  int me = brn->me;
  int nproc = brn->nproc;

  int nlocal = brn->nlocal;
  int nall = brn->nall;
  tagint nvoxel = brn->nvoxel;

  tagint *tag = brn->tag;

  dsize = num_agents + 2;
  create(send_buf,nlocal*dsize,"send_buf");

  // pack
  c = 0;
  for (i=0; i<nall; i++) {
    if(!brn->is_loc[i]) continue;
    send_buf[c++] = ubuf(tag[i]).d;
    send_buf[c++] = ubuf(brn->type[i]).d;

    for (ag_id=0; ag_id<num_agents; ag_id++)
      send_buf[c++] = brn->agent[ag_id][i];
  }

  if (!me) {
    create(rcounts,nproc,"rcounts");
    create(displs,nproc,"displs");
    create(recv_buf,nvoxel*dsize,"recv_buf");
  }

  MPI_Gather(&nlocal,1,MPI_INT,rcounts,1,MPI_INT,0,world);

  if (!me) {
    int offset = 0;
    for (i = 0; i < nproc; i++) {
      rcounts[i] *= dsize;
      displs[i] = offset;
      offset += rcounts[i];
    }
  }

  MPI_Gatherv(send_buf,nlocal*dsize,MPI_DOUBLE,recv_buf,rcounts,displs,MPI_DOUBLE,0,world);

  destroy(send_buf);

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

    c = 0;
    for (i=0; i<nvoxel; i++) {
      buft = (tagint) ubuf(recv_buf[c++]).i; // tag
      fb.write((char *)&buft,sizeof(buft));

      bufi = (int) ubuf(recv_buf[c++]).i; // type
      fb.write((char *)&bufi,sizeof(bufi));

      for (ag_id=0; ag_id<num_agents; ag_id++) {
        bufd = recv_buf[c++];
        fb.write((char *)&bufd,sizeof(bufd));
      }
    }

    fb.close();

    //read and check
    //fw = fopen("test.txt","w");
    //fb.open(dname,ios::binary|ios::in|ios::app);
    //fb.seekg(0);
    //fb.read((char *)&bufi,sizeof(bufi)); // step
    //fprintf(fw,"%i \n", bufi);
    //fb.read((char *)&buft,sizeof(buft)); //nvoxel
    //fprintf(fw,TAGINT_FORMAT " \n", buft);
    //for (i=0; i<nvoxel; i++) {
    //  fb.read((char *)&buft,sizeof(buft)); // tag
    //  fprintf(fw,TAGINT_FORMAT " ", buft);
    //  fb.read((char *)&bufi,sizeof(bufi)); // type
    //  fprintf(fw,"%i ", bufi);
    //  fb.read((char *)&bufd,sizeof(bufd)); // fAb
    //  fprintf(fw,"%g \n", bufd);
    //}
    //fb.close();
    //fclose(fw);

    destroy(rcounts);
    destroy(displs);
    destroy(recv_buf);
  }

}

/* ----------------------------------------------------------------------*/
void Output::statistics(Brain *brn) {
  if (brn->step % severy != 0) return;

  tagint i,c,nr,size;
  int ag_id;

  double *send_buf = NULL;
  double *recv_buf = NULL;

  MPI_Comm world = brn->world;

  int me = brn->me;

  double vlen = brn->vlen;
  double vlen_1 = brn->vlen_1;

  int nall = brn->nall;

  double **x = brn->x;

  double **agent_val;
  tagint *agent_num;

  nr = 2; // number of regions: parenchyma and CSF

  create(agent_val,num_agents,nr,"agent_val");
  create(agent_num,nr,"agent_num");

  for (ag_id=0; ag_id<num_agents; ag_id++)
    for (i=0; i<nr; i++) {
      agent_val[ag_id][i] = 0.0;
      agent_num[i] = 0;
    }

  for (i=0; i<nall; i++) {
    if (!brn->is_loc[i]) continue;
    if (brn->type[i] == WM_type || brn->type[i] == GM_type)
      c = 0;
    else if (brn->type[i] == CSF_type)
      c = 1;
    else
      continue;

    for (ag_id=0; ag_id<num_agents; ag_id++)
      agent_val[ag_id][c] += brn->agent[ag_id][i];
    agent_num[c]++;
  }

  size = nr*(num_agents+1);
  create(send_buf,size,"send_buf");

  // pack
  c = 0;
  for (i=0; i<nr; i++) {
    send_buf[c++] = ubuf(agent_num[i]).d;
    for (ag_id=0; ag_id<num_agents; ag_id++)
      send_buf[c++] = agent_val[ag_id][i];
  }

  if (!me)
    create(recv_buf,size,"recv_buf");

  MPI_Reduce(send_buf,recv_buf,size,MPI_DOUBLE,MPI_SUM,0,world);

  destroy(send_buf);

  if (!me) {
    // unpack
    c = 0;
    for (i=0; i<nr; i++) {
      agent_num[i] = (tagint) ubuf(recv_buf[c++]).i;
      for (ag_id=0; ag_id<num_agents; ag_id++)
        agent_val[ag_id][i] = recv_buf[c++] / agent_num[i];
    }

    // print
    fw = fopen(sname.c_str(),"a");

    // PAR
    fprintf(fw,"%i ", brn->step);
    fprintf(fw,"PAR ");
    for (ag_id=0; ag_id<num_agents; ag_id++)
      fprintf(fw,"%g ", agent_val[ag_id][0]); // agent_values
    fprintf(fw,"\n");

    // CSF
    fprintf(fw,"%i ", brn->step);
    fprintf(fw,"CSF ");
    for (ag_id=0; ag_id<num_agents; ag_id++)
      fprintf(fw,"%g ", agent_val[ag_id][1]); // agent_values
    fprintf(fw,"\n");

    fclose(fw);

    destroy(recv_buf);
  }

  destroy(agent_val);
  destroy(agent_num);

}

/* ----------------------------------------------------------------------*/
void Output::statistics_sphere(Brain *brn) {
  if (brn->step % severy != 0) return;

  tagint i,c,nr,size;
  int ag_id;

  double *send_buf = NULL;
  double *recv_buf = NULL;

  MPI_Comm world = brn->world;

  int me = brn->me;

  double vlen = brn->vlen;
  double vlen_1 = brn->vlen_1;

  int nlocal = brn->nlocal;
  int nall = brn->nall;

  double **x = brn->x;

  double **agent_val;
  tagint *agent_num;

  nr = (tagint)( brn->lbox[0] * vlen_1 / 2 );

  create(agent_val,num_agents,nr,"agent_val");
  create(agent_num,nr,"agent_num");

  for (ag_id=0; ag_id<num_agents; ag_id++)
    for (i=0; i<nr; i++) {
      agent_val[ag_id][i] = 0.0;
      agent_num[i] = 0;
    }

  for (i=0; i<nall; i++) {
    if(!brn->is_loc[i]) continue;
    c = (tagint) (vlen_1 * sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1] + x[i][2] * x[i][2]));
    //c = (tagint) (vlen_1 * (x[i][0] - brn->boxlo[0]));
    if (c >= nr) continue;

    for (ag_id=0; ag_id<num_agents; ag_id++)
      agent_val[ag_id][c] += brn->agent[ag_id][i];
    agent_num[c]++;
  }

  size = nr*(num_agents+1);
  create(send_buf,size,"send_buf");

  // pack
  c = 0;
  for (i=0; i<nr; i++) {
    send_buf[c++] = ubuf(agent_num[i]).d;
    for (ag_id=0; ag_id<num_agents; ag_id++)
      send_buf[c++] = agent_val[ag_id][i];
  }

  if (!me)
    create(recv_buf,size,"recv_buf");

  MPI_Reduce(send_buf,recv_buf,size,MPI_DOUBLE,MPI_SUM,0,world);

  destroy(send_buf);

  if (!me) {
    // unpack
    c = 0;
    for (i=0; i<nr; i++) {
      agent_num[i] = (tagint) ubuf(recv_buf[c++]).i;
      for (ag_id=0; ag_id<num_agents; ag_id++)
        agent_val[ag_id][i] = recv_buf[c++] / agent_num[i];
    }

    // print
    fw = fopen(sname.c_str(),"a");

    fprintf(fw,"ITEM: TIMESTEP \n");
    fprintf(fw,"%i \n", brn->step);

    fprintf(fw,"ITEM: NUMBER OF BINS \n");
    fprintf(fw,TAGINT_FORMAT " \n", nr);

    fprintf(fw,"ITEM: BIN id r ");
    for (ag_id=0; ag_id<num_agents; ag_id++)
      fprintf(fw,"%s ", ag_str[ag_id].c_str());
    fprintf(fw,"\n");

    for (i=0; i<nr; i++) {
      fprintf(fw,TAGINT_FORMAT " ", i); // id
      fprintf(fw,"%g ", vlen * i); // r
      for (ag_id=0; ag_id<num_agents; ag_id++)
        fprintf(fw,"%g ", agent_val[ag_id][i]); // agent_values
      fprintf(fw,"\n");
    }

    fclose(fw);

    destroy(recv_buf);
  }

  destroy(agent_val);
  destroy(agent_num);

}

/* ----------------------------------------------------------------------*/
void Output::dump(Brain *brn) {
  int i;

  for (i=0; i<dump_arg.size(); i++) {

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
void Output::dump_txt(Brain *brn, vector<string> arg) {
  tagint i,c;
  int dsize, ag_id, aid;

  int *rcounts = NULL;
  int *displs = NULL;
  double *send_buf = NULL;
  double *recv_buf = NULL;

  MPI_Comm world = brn->world;

  int me = brn->me;
  int nproc = brn->nproc;

  int nlocal = brn->nlocal;
  int nall = brn->nall;
  tagint nvoxel = brn->nvoxel;

  double **x = brn->x;

  dsize = 3;
  c = 3;
  while (c < arg.size()) {
    if (!arg[c].compare("type"))
      dsize++;
    else if (brn->input->find_agent(arg[c]) >= 0)
      dsize++;
    else {
      printf("Error: dump argument \"%s\" unknown. \n", arg[c].c_str());
      exit(1);
    }
    c++;
  }

  create(send_buf,nlocal*dsize,"send_buf");

  // pack
  c = 0;
  for (i=0; i<nall; i++) {
    if (!brn->is_loc[i]) continue;
    send_buf[c++] = x[i][0]; // x
    send_buf[c++] = x[i][1]; // y
    send_buf[c++] = x[i][2]; // z

    aid = 3;
    while (aid < arg.size()) {
      ag_id = brn->input->find_agent(arg[aid]);
      if (!arg[aid].compare("type"))
        send_buf[c++] = ubuf(brn->type[i]).d;
      else if (ag_id >= 0)
        send_buf[c++] = brn->agent[ag_id][i];
      aid++;
    }

  }

  if (!me) {
    create(rcounts,nproc,"rcounts");
    create(displs,nproc,"displs");
    create(recv_buf,nvoxel*dsize,"recv_buf");
  }

  MPI_Gather(&nlocal,1,MPI_INT,rcounts,1,MPI_INT,0,world);

  if (!me) {
    int offset = 0;
    for (i = 0; i < nproc; i++) {
      rcounts[i] *= dsize;
      displs[i] = offset;
      offset += rcounts[i];
    }
  }

  MPI_Gatherv(send_buf,nlocal*dsize,MPI_DOUBLE,recv_buf,rcounts,displs,MPI_DOUBLE,0,world);

  destroy(send_buf);

  // unpack and print on the root
  if (!me) {
    sort_tag(brn,recv_buf,dsize);

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

    aid = 3;
    while (aid < arg.size()) {
      ag_id = brn->input->find_agent(arg[aid]);
      if (!arg[aid].compare("type"))
        fprintf(fw,"type ");
      else if (ag_id >= 0)
        fprintf(fw,"%s ", ag_str[ag_id].c_str());
      aid++;
    }
    fprintf(fw,"\n");

    c = 0;
    for (i=0; i<nvoxel; i++) {
      c += 3; // skip x,y,z
      fprintf(fw,TAGINT_FORMAT " ", i); // tag

      aid = 3;
      while (aid < arg.size()) {
        ag_id = brn->input->find_agent(arg[aid]);
        if (!arg[aid].compare("type"))
          fprintf(fw,"%i ", (int) ubuf(recv_buf[c++]).i); // type
        else if (ag_id >= 0)
          fprintf(fw,"%g ", recv_buf[c++]); // agent
        aid++;
      }
      fprintf(fw,"\n");

    }

    fclose(fw);

    destroy(rcounts);
    destroy(displs);
    destroy(recv_buf);
  }

}

/* ----------------------------------------------------------------------*/
void Output::dump_mri(Brain *brn, vector<string> arg) {
  tagint i,c;
  int dsize, ag_id, aid;

  int *rcounts = NULL;
  int *displs = NULL;
  double *send_buf = NULL;
  double *recv_buf = NULL;

  MPI_Comm world = brn->world;

  int me = brn->me;
  int nproc = brn->nproc;

  int nlocal = brn->nlocal;
  int nall = brn->nall;
  tagint nvoxel = brn->nvoxel;

  double **x = brn->x;

  dsize = 3;
  c = 3;
  while (c < arg.size()) {
    if (!arg[c].compare("type"))
      dsize++;
    else if (!arg[c].compare("me"))
      dsize++;
    else if (brn->input->find_agent(arg[c]) >= 0)
      dsize++;
    else {
      printf("Error: dump argument \"%s\" unknown. \n", arg[c].c_str());
      exit(1);
    }
    c++;
  }

  create(send_buf,nlocal*dsize,"send_buf");

  // pack
  c = 0;
  for (i=0; i<nall; i++) {
    if(!brn->is_loc[i]) continue;
    send_buf[c++] = x[i][0]; // x
    send_buf[c++] = x[i][1]; // y
    send_buf[c++] = x[i][2]; // z

    aid = 3;
    while (aid < arg.size()) {
      ag_id = brn->input->find_agent(arg[aid]);
      if (!arg[aid].compare("type"))
        send_buf[c++] = ubuf(brn->type[i]).d;
      else if (!arg[aid].compare("me"))
        send_buf[c++] = ubuf(me).d;
      else if (ag_id >= 0)
        send_buf[c++] = brn->agent[ag_id][i];
      aid++;
    }
  }

  if (!me) {
    create(rcounts,nproc,"rcounts");
    create(displs,nproc,"displs");
    create(recv_buf,nvoxel*dsize,"recv_buf");
  }

  MPI_Gather(&nlocal,1,MPI_INT,rcounts,1,MPI_INT,0,world);

  if (!me) {
    int offset = 0;
    for (i = 0; i < nproc; i++) {
      rcounts[i] *= dsize;
      displs[i] = offset;
      offset += rcounts[i];
    }
  }

  MPI_Gatherv(send_buf,nlocal*dsize,MPI_DOUBLE,recv_buf,rcounts,displs,MPI_DOUBLE,0,world);

  destroy(send_buf);

  // unpack and print on the root
  if (!me) {
    sort_tag(brn,recv_buf,dsize);

    const int dims5[] = {5, brn->nv[0], brn->nv[1], brn->nv[2], 1, dsize-3, 1, 1};
    const int dims3[] = {3, brn->nv[0], brn->nv[1], brn->nv[2], 1, 1, 1, 1};

    nifti_image *nim;

    if (dsize > 4)
      nim = nifti_image_setup(brn,arg, dims5, NIFTI_INTENT_VECTOR);
    else
      nim = nifti_image_setup(brn,arg, dims3, NIFTI_INTENT_NONE);

    float* data = (float*) nim->data;

    tagint cnim = 0;
    c = 0;
    for (int kk=0; kk<nim->nz; kk++)
      for (int jj=0; jj<nim->ny; jj++)
        for (int ii=0; ii<nim->nx; ii++) {
          c += 3;

          aid = 3;
          while (aid < arg.size()) {
            cnim = ii + nim->nx * ( jj + nim->ny * (kk + nim->nz * ( (aid-3) ) ) );
            //cnim = (aid-3) + (arg.size()-3) * ( kk + nim->nz * (jj + nim->ny * ii ) );
            ag_id = brn->input->find_agent(arg[aid]);

            if (!arg[aid].compare("type"))
              data[cnim] = (float) ubuf(recv_buf[c++]).i;
            else if (!arg[aid].compare("me"))
              data[cnim] = (float) ubuf(recv_buf[c++]).i;
            else if (ag_id >= 0)
              data[cnim] = (float) recv_buf[c++];

            //if (!arg[aid].compare("me"))
             // printf("proc %i: HERE0 output itag=%i, data[%i, %i, %i, %i] = %g \n",
               //      brn->me, cnim, ii,jj,kk,aid-3, data[cnim]);
            //printf("proc %i: HERE1 output itag=%i, data[%i, %i, %i, %i] = %g \n",
              //     brn->me, itag, i,j,k,h, (float) ptr[c]);

            aid++;
          }
        }

    nifti_image_write(nim);
    nifti_image_free(nim);

    destroy(rcounts);
    destroy(displs);
    destroy(recv_buf);
  }

}

/* ----------------------------------------------------------------------*/
void Output::sort_tag(Brain *brn, double *data, int dsize) {
  tagint i, ii, itag, it;
  double xx[3];
  double dum;

  i = 0;
  while (i < brn->nvoxel) {
    ii = i*dsize;
    xx[0] = data[ii];
    xx[1] = data[ii+1];
    xx[2] = data[ii+2];

  //printf("proc %i: HERE2, %i %i %i %g %g %g \n",brn->me,i,ii,dsize,brn->boxlo[0],brn->boxlo[1],brn->boxlo[2]);

    itag = find_tag(brn, xx[0], xx[1], xx[2]);

    if (i == itag)
      i++;
    else {
      // exchange data
      it = itag*dsize;
      for (int j=0; j<dsize; j++) {
        dum = data[it + j];
        data[it + j] = data[ii + j];
        data[ii + j] = dum;
      }
    }

  }

}

/* ----------------------------------------------------------------------*/
tagint Output::find_tag(Brain *brn, double x, double y, double z) {
  int i,j,k;
  tagint itag;

  int *nv = brn->nv;
  double *boxlo = brn->boxlo;
  double vlen_1 = brn->vlen_1;

  i = (int) ((x - boxlo[0]) * vlen_1 - 0.5);
  j = (int) ((y - boxlo[1]) * vlen_1 - 0.5);
  k = (int) ((z - boxlo[2]) * vlen_1 - 0.5);

  itag = i + nv[0] * (j + nv[1]*k);

//printf("proc %i: HERE1 "TAGINT_FORMAT", %g %g %g %g %g %g \n",brn->me,itag,x,y,z, boxlo[0],boxlo[1],boxlo[2]);

  return itag;
}

/* ----------------------------------------------------------------------*/
nifti_image* Output::nifti_image_setup(Brain *brn, vector<string> arg,
                                       const int dims[], int intent) {
  int i;
  string str;

  nifti_image *nim;
  //nifti_1_header *nhdr;

  //create(dims,8,"output:dims");

  //nhdr = nifti_make_new_header(dims, DT_FLOAT32);
  nim = nifti_make_new_nim(dims,DT_FLOAT32,1);

  // header and image file names (.nii)
  str = "";
  str.append(arg[2]);
  str.append(to_string(brn->step));
  const char *prefix = str.c_str();

  if( nifti_set_filenames(nim, prefix, 1, 1) ) {
    printf("Error: nifti_set_filenames cannot set the names for %s . \n", prefix);
    exit(1);
  }

  //nhdr->sizeof_hdr = 348;

  for (i=0; i<8; i++)
    nim->dim[i] = dims[i];

  nim->intent_p1 = nim->intent_p2 = nim->intent_p3 = 0.0;

  nim->intent_code = intent; // NIFTI_INTENT_VECTOR;

  nim->datatype = DT_FLOAT32;
  nim->nbyper = 16;

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
  for (i=1; i<8; i++)
    nim->nvox *= nim->dim[i];

  //nhdr->vox_offset = 352;

  nim->scl_slope = 0.0;
  nim->scl_inter = 0.0;

  nim->slice_end = 0.0;
  nim->slice_code = 0.0;

  nim->cal_min = 0.0;
  nim->cal_max = 0.0;

  //nhdr->xyzt_units = NIFTI_UNITS_MICRON;
  nim->xyz_units = NIFTI_UNITS_MICRON;
  nim->time_units = NIFTI_UNITS_UNKNOWN;

  nim->slice_duration = 0.0;
  // initial time of the image
  nim->toffset = brn->step * brn->dt;

  //nhdr->glmax = INT_MAX; // unused
  //nhdr->glmin = INT_MIN; // unused

  str = "";
  for (i=3; i<arg.size(); i++) {
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

/*
/* ----------------------------------------------------------------------
nifti_image* Output::nifti_image_setup(Brain *brn, vector<string> arg, int dsize) {
  int i;

  int *nv = brn->nv;

  nifti_image *nim = brn->nim;

  nifti_image* nifti_convert_nhdr2nim(struct nifti_1_header nhdr,
                                      const char * fname);

  nim->ndim = 5;
  nim->nx = nv[0];
  nim->ny = nv[1];
  nim->nz = nv[2];
  nim->nt = 1;
  nim->nu = dsize;
  nim->nv = 1;
  nim->nw = 1;

  nim->dim[0] = nim->ndim;
  nim->dim[1] = nim->nx;
  nim->dim[2] = nim->ny;
  nim->dim[3] = nim->nz;
  nim->dim[4] = nim->nt;
  nim->dim[5] = nim->nu;
  nim->dim[6] = nim->nv;
  nim->dim[7] = nim->nw;

  nim->nvox = 1;
  for (i=1; i<8; i++)
    nim->nvox *= nim->dim[i];

  nim->nbyper = 16;
  nim->datatype = DT_FLOAT;

  nim->dx = brn->vlen;
  nim->dy = brn->vlen;
  nim->dz = brn->vlen;
  nim->dt = brn->dt * stoi(arg[1]);
  nim->du = 1;
  nim->dv = 1;
  nim->dw = 1;

  nim->pixdim[0] = 1;
  nim->pixdim[1] = nim->dx;
  nim->pixdim[2] = nim->dy;
  nim->pixdim[3] = nim->dz;
  nim->pixdim[4] = nim->dt;
  nim->pixdim[5] = nim->du;
  nim->pixdim[6] = nim->dv;
  nim->pixdim[7] = nim->dw;

  nim->scl_slope = 0.0;
  nim->scl_inter = 0.0;

  nim->cal_min = 0.0;
  nim->cal_max = 0.0;

  // voxel spatial position mapping scheme
  nim->qform_code = 0; // Method 2: using quaternions
  nim->sform_code = 0; // Method 3: using general affine transformation

  // not so important parameters
  nim->freq_dim = 0;
  nim->phase_dim = 0;
  nim->slice_dim = 0;
  nim->slice_code = 0;
  nim->slice_start = 0;
  nim->slice_end = 0;
  nim->slice_duration = 0.0;

  // quaternions - Method 2 of mapping
  nim->quatern_b = nim->quatern_c = nim->quatern_d = 0.0;
  nim->qoffset_x = nim->qoffset_y = nim->qoffset_z = 0.0;
  nim->qfac = 1.0;

  // initial time of the image
  nim->toffset = brn->step * brn->dt;

  nim->xyz_units = NIFTI_UNITS_MICRON;
  nim->time_units = NIFTI_UNITS_UNKNOWN;

  //nim->nifti_type = 1; // NIFTI-1 (1 file)

  nim->intent_code = NIFTI_INTENT_VECTOR;
  nim->intent_p1 = nim->intent_p2 = nim->intent_p3 = 0.0;

  string str = "";
  for (i=3; i<arg.size(); i++)
    str.append(arg[i]);

  strcpy(nim->descrip, str.c_str());

  memset(nim->intent_name, 0, sizeof(nim->intent_name));
  memset(nim->aux_file, 0, sizeof(nim->aux_file));

  // header and image file names (.nii)
  //str = "";
  //str.append(arg[2]);
  //str.append(to_string(brn->step));
  //str.append(".nii");

  //nim->fname = (char*) str.c_str();
  //nim->iname = (char*) str.c_str();

  //  int swapsize ;              /*!< swap unit in image data (might be 0)
  //  int byteorder ;             /*!< byte order on disk (MSB_ or LSB_FIRST)
  if (nim->data != NULL) {
    destroy(nim->data);
    nim->data = NULL;
  }

  // allocate memory for nifti data
  nim->data = smalloc(nim->nbyper*nim->nvox,"output:nifti_image_data");

  nim->num_ext = 0;

  return nim;

}
*/
