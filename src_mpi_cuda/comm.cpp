#include <mpi.h>
#include "math.h"

#include "pointers.h"
#include "brain.h"
#include "comm.h"

using namespace std;
using namespace brain_NS;

/* ---------------------------------------------------------------------- */
Comm::Comm() {
  b_itr = 0;
  b_dim.assign("none");

  comm_size = 2 + num_agents;
}

/* ----------------------------------------------------------------------*/
Comm::~Comm() {
  destroy(send_buf);
  destroy(recv_buf);
}

/* ----------------------------------------------------------------------
 * Set partitions for the parallel simulation.
 * For each partition: find boundaries xlo[i] and xhi[i],
 * and its communicating neighbor partitions.
 * ----------------------------------------------------------------------*/
void Comm::partition(Brain *brn) {
  int i,j,k;
  double pos[3],POS[3],DX[3];

  int *npart = brn->npart;

  double *boxlo = brn->boxlo;

  double *xlo = brn->xlo;
  double *xhi = brn->xhi;

  double *lbox = brn->lbox;

  for (i=0; i<3; i++)
    DX[i] = lbox[i] / npart[i];

  int me = 0;
  for (i=0; i<npart[0]; i++) {
    pos[0] = boxlo[0] + 0.5*DX[0] + DX[0]*i;
    for (j=0; j<npart[1]; j++) {
      pos[1] = boxlo[1] + 0.5*DX[1] + DX[1]*j;
      for (k=0; k<npart[2]; k++) {
        pos[2] = boxlo[2] + 0.5*DX[2] + DX[2]*k;
        if (me == brn->me) {
          POS[0] = pos[0];
          POS[1] = pos[1];
          POS[2] = pos[2];

          // set comm_side
          comm_side[XLO] = find_me(brn,i-1,j,k);
          comm_side[XHI] = find_me(brn,i+1,j,k);
          comm_side[YLO] = find_me(brn,i,j-1,k);
          comm_side[YHI] = find_me(brn,i,j+1,k);
          comm_side[ZLO] = find_me(brn,i,j,k-1);
          comm_side[ZHI] = find_me(brn,i,j,k+1);
        }
        me++;
      }
    }
  }

  for (i=0; i<3; i++) {
    xlo[i] = POS[i] - 0.5*DX[i];
    xhi[i] = POS[i] + 0.5*DX[i];
  }

}

/* ----------------------------------------------------------------------
 * Balance the partition loop sizes among all the partitions in a specified
 * dimension, dim;
 * ----------------------------------------------------------------------*/
void Comm::balance(Brain *brn) {
  int b_flag, cid, sid;

  if (!b_dim.compare("x"))
    b_flag = 0;
  else if (!b_dim.compare("y"))
    b_flag = 1;
  else if ((!b_dim.compare("z")))
    b_flag = 2;
  else {
    printf("Error: balance dimension unknown. \n");
    exit(1);
  }

  int i,j,k;
  int nloop_loc;

  int me = brn->me;
  int nproc = brn->nproc;

  MPI_Comm world = brn->world;

  int *nv = brn->nv;

  int *npart = brn->npart;
  double *xlo = brn->xlo;
  double *xhi = brn->xhi;

  int nall = brn->nall;

  int *type = brn->type;

  nloop_loc = 0;
  for (i=0; i<nall; i++) {
    if (type[i] == EMP_type) continue;
    nloop_loc++;
  }

  int dsize = 3; // data size for communication; nloop, xlo, and xhi

  int *rcounts = NULL;
  int *displs = NULL;
  double *r_buf = NULL;

  create(rcounts,nproc,"comm:rcounts");
  create(displs,nproc,"comm:displs");
  create(r_buf,nproc*dsize,"comm:r_buf");

  double *s_buf = NULL;
  create(s_buf,dsize,"comm:s_buf");

  // pack
  s_buf[0] = ubuf(nloop_loc).d;
  s_buf[1] = xlo[b_flag];
  s_buf[2] = xhi[b_flag];

  int offset = 0;
  for (i = 0; i < nproc; i++) {
    rcounts[i] = dsize;
    displs[i] = offset;
    offset += rcounts[i];
  }

  MPI_Allgatherv(s_buf,dsize,MPI_DOUBLE,r_buf,rcounts,displs,MPI_DOUBLE,world);

  destroy(rcounts);
  destroy(displs);

  destroy(s_buf);

  // unpack
  int *nloop = NULL;
  double *xlo_all = NULL;
  double *xhi_all = NULL;

  create(nloop,nproc,"comm:nloop");
  create(xlo_all,nproc,"comm:xlo_all");
  create(xhi_all,nproc,"comm:xhi_all");

  int c = 0;
  for (i=0; i<nproc; i++) {
    c = i*dsize;
    nloop[i] = (int) ubuf(r_buf[c]).i;
    xlo_all[i] = r_buf[c+1];
    xhi_all[i] = r_buf[c+2];
  }

  destroy(r_buf);

  double *nloop_sec = NULL;
  double *xlo_sec = NULL;
  double *xhi_sec = NULL;
  int *n_sec = NULL;

  create(nloop_sec,npart[b_flag],"comm:nloop_sec");
  create(xlo_sec,npart[b_flag],"comm:xlo_sec");
  create(xhi_sec,npart[b_flag],"comm:xhi_sec");
  create(n_sec,npart[b_flag],"comm:n_sec");

  for (i=0; i<npart[b_flag]; i++) {
    n_sec[i] = 0;
    nloop_sec[i] = 0.0;
  }

  for (i=0; i<npart[0]; i++) {
    if (b_flag == 0)
      sid = i;
    for (j=0; j<npart[1]; j++) {
      if (b_flag == 1)
        sid = j;
      for (k=0; k<npart[2]; k++) {
        if (b_flag == 2)
          sid = k;

        cid = find_me(brn,i,j,k);
        nloop_sec[sid] += (double)(nloop[cid]);
        xlo_sec[sid] = xlo_all[cid];
        xhi_sec[sid] = xhi_all[cid];
        n_sec[sid]++;
      }
    }
  }

  double nloop_avg = 0.0;
  for (i=0; i<npart[b_flag]; i++) {
    nloop_sec[i] /= n_sec[i];
    nloop_avg += nloop_sec[i];
  }

  nloop_avg /= npart[b_flag];

  double lx, lx_new, correction;

  double *xlo_new = NULL;
  double *xhi_new = NULL;

  create(xlo_new,npart[b_flag],"comm:xlo_new");
  create(xhi_new,npart[b_flag],"comm:xhi_new");

  double dum = (double)(nv[b_flag]) / brn->nvoxel * brn->vlen;

  xlo_new[0] = xlo_sec[0];
  for (i=0; i<npart[b_flag]; i++) {
    correction = (nloop_avg - nloop_sec[i]) * dum;
    lx = xhi_sec[i] - xlo_sec[i];
    lx_new = lx + correction;
    if (i > 0)
      xlo_new[i] = xhi_new[i-1];
    xhi_new[i] = xlo_new[i] + lx_new;
  }
  xhi_new[i] = xhi_sec[i];

  // set the new xlo and xhi
  int dumi;
  i = (int) (me / (npart[1]*npart[2]));
  dumi = me % (npart[1]*npart[2]);
  j = (int) (dumi / npart[2]);
  k = dumi % npart[2];

  if (b_flag == 0)
    sid = i;
  if (b_flag == 1)
    sid = j;
  if (b_flag == 2)
    sid = k;

  //for (i=0; i<npart[b_flag]; i++)
  //  printf("proc %i: HERE sec %i, %g %g  %g %g \n",me,i,xlo[b_flag],xhi[b_flag],xlo_new[i],xhi_new[i]);

  xlo[b_flag] = xlo_new[sid];
  xhi[b_flag] = xhi_new[sid];

  // write result
  int nloop_max, nloop_min;

  //printf("proc %i: HERE1 %i %i \n",me,nloop_loc,nloop_max);
  MPI_Reduce(&nloop_loc,&nloop_max,1,MPI_INT,MPI_MAX,0,world);
  MPI_Reduce(&nloop_loc,&nloop_min,1,MPI_INT,MPI_MIN,0,world);

  if (!me)
    printf("Balance: unevenness ratio (max/min) before balancing = %g \n",
           (float)(nloop_max)/nloop_min);

  //printf("proc %i: HERE %g %g  %g %g  %g %g \n",me, xlo[0], xhi[0], xlo[1], xhi[1], xlo[2], xhi[2]);

  destroy(nloop);
  destroy(xlo_all);
  destroy(xhi_all);

  destroy(nloop_sec);
  destroy(xlo_sec);
  destroy(xhi_sec);
  destroy(n_sec);

  destroy(xlo_new);
  destroy(xhi_new);

}

/* ----------------------------------------------------------------------
 * Initialize communications: each partition has 6 communicating neighbors,
 * 3 for forward communication and 3 for reverse communication.
 * set maximum buffer size for communications, and the bufer size for each
 * communication.
 * ----------------------------------------------------------------------*/
void Comm::comm_init(Brain *brn) {
  int i,flag;

  int nall = brn->nall;

  int *is_loc = brn->is_loc;
  double **x = brn->x;

  double *xlo = brn->xlo;
  double *xhi = brn->xhi;

  // find buffer size
  max_buf_size = 0;
  for (flag=0; flag<6; flag++) {
    int c = 0;
    for (i=0; i<nall; i++) {
      if (is_loc[i]) continue;
      if (flag == XLO && x[i][0] < xlo[0]) c++;
      else if (flag == XHI && x[i][0] >= xhi[0]) c++;
      else if (flag == YLO && x[i][1] < xlo[1]) c++;
      else if (flag == YHI && x[i][1] >= xhi[1]) c++;
      else if (flag == ZLO && x[i][2] < xlo[2]) c++;
      else if (flag == ZHI && x[i][2] >= xhi[2]) c++;
    }

    if (c > max_buf_size)
      max_buf_size = c;
    buf_size[flag] = c * comm_size;
  }

  max_buf_size *= comm_size;

  int dum = max_buf_size + 1; // the first element of buffer is buf_size
  // set max of max_buf_size as a universal communication size
  MPI_Allreduce(&dum,&max_buf_size,1,MPI_INT,MPI_MAX,brn->world);

  allocations(brn);

}

/* ----------------------------------------------------------------------*/
void Comm::allocations(Brain *brn) {
  create(send_buf,max_buf_size,"send_buf");
  create(recv_buf,max_buf_size,"recv_buf");
}

/* ----------------------------------------------------------------------
 * Find the rank of a partition, with coordinates i,j,k
 * ----------------------------------------------------------------------*/
int Comm::find_me(Brain *brn, int i, int j, int k) {
  int *npart = brn->npart;

  if (i < 0 || i >= npart[0])
    return -1;
  if (j < 0 || j >= npart[1])
    return -1;
  if (k < 0 || k >= npart[2])
    return -1;

  return k + npart[2] * (j + npart[1]*i);

}

/* ----------------------------------------------------------------------
 * Forward communication: communicate the properties from the local voxels
 * to the ghost voxels of neighboring partitions.
 * ----------------------------------------------------------------------*/
void Comm::forward_comm(Brain *brn) {
  int ctag;
  MPI_Request req_send,req_recv;
  MPI_Comm world = brn->world;

  ctag = 0;

  /// X direction
  if (comm_side[XHI] >= 0)
    MPI_Irecv(recv_buf,max_buf_size,MPI_DOUBLE,comm_side[XHI],ctag,world,&req_recv);

  if (comm_side[XLO] >= 0) {
    forward_pack(brn,XLO);
    MPI_Isend(send_buf,max_buf_size,MPI_DOUBLE,comm_side[XLO],ctag,world,&req_send);
  }

  if (comm_side[XHI] >= 0) {
    MPI_Wait(&req_recv,MPI_STATUS_IGNORE);
    forward_unpack(brn);
  }

  if (comm_side[XLO] >= 0)
    MPI_Wait(&req_send,MPI_STATUS_IGNORE);

  /// Y direction
  if (comm_side[YHI] >= 0)
    MPI_Irecv(recv_buf,max_buf_size,MPI_DOUBLE,comm_side[YHI],ctag,world,&req_recv);

  if (comm_side[YLO] >= 0) {
    forward_pack(brn,YLO);
    MPI_Isend(send_buf,max_buf_size,MPI_DOUBLE,comm_side[YLO],ctag,world,&req_send);
  }

  if (comm_side[YHI] >= 0) {
    MPI_Wait(&req_recv,MPI_STATUS_IGNORE);
    forward_unpack(brn);
  }

  if (comm_side[YLO] >= 0)
    MPI_Wait(&req_send,MPI_STATUS_IGNORE);


  /// Z direction
  if (comm_side[ZHI] >= 0)
    MPI_Irecv(recv_buf,max_buf_size,MPI_DOUBLE,comm_side[ZHI],ctag,world,&req_recv);

  if (comm_side[ZLO] >= 0) {
    forward_pack(brn,ZLO);
    MPI_Isend(send_buf,max_buf_size,MPI_DOUBLE,comm_side[ZLO],ctag,world,&req_send);
  }

  if (comm_side[ZHI] >= 0) {
    MPI_Wait(&req_recv,MPI_STATUS_IGNORE);
    forward_unpack(brn);
  }

  if (comm_side[ZLO] >= 0)
    MPI_Wait(&req_send,MPI_STATUS_IGNORE);

}

/* ----------------------------------------------------------------------
 * Packing the buffer for forward communication
 * ----------------------------------------------------------------------*/
void Comm::forward_pack(Brain *brn, int flag) {
  int i, ag_id;

  double vlen = brn->vlen;
  int nlocal = brn->nlocal;

  int *nvl = brn->nvl;

  //double *xlo = brn->xlo;
  //double **x = brn->x;

  int c = 0;
  send_buf[c++] = ubuf(buf_size[flag]).d;

  if (flag == XLO) {
    int ii = 0;
    for (int jj=0; jj<nvl[1]; jj++)
      for (int kk=0; kk<nvl[2]; kk++) {
        i = kk + nvl[2] * (jj + nvl[1]*ii);

        send_buf[c++] = ubuf(brn->tag[i]).d;
        send_buf[c++] = ubuf(brn->type[i]).d;
        for (ag_id=0; ag_id<num_agents; ag_id++)
          send_buf[c++] = brn->agent[ag_id][i];
      }
  }

  else if (flag == YLO) {
    int jj = 0;
    for (int ii=0; ii<nvl[0]; ii++)
      for (int kk=0; kk<nvl[2]; kk++) {
        i = kk + nvl[2] * (jj + nvl[1]*ii);
   
        send_buf[c++] = ubuf(brn->tag[i]).d;
        send_buf[c++] = ubuf(brn->type[i]).d;
        for (ag_id=0; ag_id<num_agents; ag_id++)
          send_buf[c++] = brn->agent[ag_id][i];
      }
  }

  else if (flag == ZLO) {
    int kk = 0;
    for (int ii=0; ii<nvl[0]; ii++)
      for (int jj=0; jj<nvl[1]; jj++) {
        i = kk + nvl[2] * (jj + nvl[1]*ii);
    
        send_buf[c++] = ubuf(brn->tag[i]).d;
        send_buf[c++] = ubuf(brn->type[i]).d;
        for (ag_id=0; ag_id<num_agents; ag_id++)
          send_buf[c++] = brn->agent[ag_id][i];
      }
  }

  else if (flag == XHI) {
    int ii = nvl[0] - 1;
    for (int jj=0; jj<nvl[1]; jj++)
      for (int kk=0; kk<nvl[2]; kk++) {
        i = kk + nvl[2] * (jj + nvl[1]*ii);

        send_buf[c++] = ubuf(brn->tag[i]).d;
        send_buf[c++] = ubuf(brn->type[i]).d;
        for (ag_id=0; ag_id<num_agents; ag_id++)
          send_buf[c++] = brn->agent[ag_id][i];
      }
  }

  else if (flag == YHI) {
    int jj = nvl[1] - 1;
    for (int ii=0; ii<nvl[0]; ii++)
      for (int kk=0; kk<nvl[2]; kk++) {
        i = kk + nvl[2] * (jj + nvl[1]*ii);

        send_buf[c++] = ubuf(brn->tag[i]).d;
        send_buf[c++] = ubuf(brn->type[i]).d;
        for (ag_id=0; ag_id<num_agents; ag_id++)
          send_buf[c++] = brn->agent[ag_id][i];
      }
  }

  else if (flag == ZHI) {
    int kk = nvl[2] - 1;
    for (int ii=0; ii<nvl[0]; ii++)
      for (int jj=0; jj<nvl[1]; jj++) {
        i = kk + nvl[2] * (jj + nvl[1]*ii);

        send_buf[c++] = ubuf(brn->tag[i]).d;
        send_buf[c++] = ubuf(brn->type[i]).d;
        for (ag_id=0; ag_id<num_agents; ag_id++)
          send_buf[c++] = brn->agent[ag_id][i];
      }
  }

  //for (i=0; i<nlocal; i++) {
    //if (flag == XLO && x[i][0] >= xlo[0] + vlen) continue;
    //else if (flag == YLO && x[i][1] >= xlo[1] + vlen) continue;
    //else if (flag == ZLO && x[i][2] >= xlo[2] + vlen) continue;

    //printf("proc %i: pack itag = %li, c=%i, i=%i, buf_size =%i, flag=%i \n",
          // brn->me,brn->tag[i],c,i,buf_size[flag],flag);

    //send_buf[c++] = ubuf(brn->tag[i]).d;
    //send_buf[c++] = ubuf(brn->type[i]).d;
    //for (ag_id=0; ag_id<num_agents; ag_id++)
      //send_buf[c++] = brn->agent[ag_id][i][0];
  //}

}

/* ----------------------------------------------------------------------
 * Unpacking the buffer for forward communication
 * ----------------------------------------------------------------------*/
void Comm::forward_unpack(Brain *brn) {
  int i, ag_id;
  tagint itag;

  int c = 0;
  int size_unpack = (int) ubuf(recv_buf[c++]).i;

  while (c <= size_unpack) {
    itag = (tagint) ubuf(recv_buf[c++]).i;
    i = brn->map[itag];

    // //////DEBUG/////////////////////
    //if (brn->me == 0) {
    //  printf("HERE %i " TAGINT_FORMAT " size_unpack=%i \n",i,itag, size_unpack);
    //}
    // //////DEBUG/////////////////////

    //printf("proc %i: unpack itag = %li, c=%i, i=%i, buf_size =%i \n",
      //     brn->me,itag,c,i,size_unpack);

    if (i != -1) {
      brn->type[i] = (int) ubuf(recv_buf[c++]).i;
      for (ag_id=0; ag_id<num_agents; ag_id++)
        brn->agent[ag_id][i] = recv_buf[c++];
    } else
      c += num_agents + 1;
  }

}

/* ----------------------------------------------------------------------
 * Reverse communication: communicate the properties from the ghost voxels
 * to the local voxels of neighboring partitions.
 * ----------------------------------------------------------------------*/
void Comm::reverse_comm(Brain *brn) {
  int ctag;
  MPI_Request req_send,req_recv;
  MPI_Comm world = brn->world;

  ctag = 0;

  /// X direction
  if (comm_side[XLO] >= 0)
    MPI_Irecv(recv_buf,max_buf_size,MPI_DOUBLE,comm_side[XLO],ctag,world,&req_recv);

  if (comm_side[XHI] >= 0) {
    reverse_pack(brn,XHI);
    MPI_Isend(send_buf,max_buf_size,MPI_DOUBLE,comm_side[XHI],ctag,world,&req_send);
  }

  if (comm_side[XLO] >= 0) {
    MPI_Wait(&req_recv,MPI_STATUS_IGNORE);
    reverse_unpack(brn);
  }

  if (comm_side[XHI] >= 0)
    MPI_Wait(&req_send,MPI_STATUS_IGNORE);

  /// Y direction
  if (comm_side[YLO] >= 0)
    MPI_Irecv(recv_buf,max_buf_size,MPI_DOUBLE,comm_side[YLO],ctag,world,&req_recv);

  if (comm_side[YHI] >= 0) {
    reverse_pack(brn,YHI);
    MPI_Isend(send_buf,max_buf_size,MPI_DOUBLE,comm_side[YHI],ctag,world,&req_send);
  }

  if (comm_side[YLO] >= 0) {
    MPI_Wait(&req_recv,MPI_STATUS_IGNORE);
    reverse_unpack(brn);
  }

  if (comm_side[YHI] >= 0)
    MPI_Wait(&req_send,MPI_STATUS_IGNORE);


  /// Z direction
  if (comm_side[ZLO] >= 0)
    MPI_Irecv(recv_buf,max_buf_size,MPI_DOUBLE,comm_side[ZLO],ctag,world,&req_recv);

  if (comm_side[ZHI] >= 0) {
    reverse_pack(brn,ZHI);
    MPI_Isend(send_buf,max_buf_size,MPI_DOUBLE,comm_side[ZHI],ctag,world,&req_send);
  }

  if (comm_side[ZLO] >= 0) {
    MPI_Wait(&req_recv,MPI_STATUS_IGNORE);
    reverse_unpack(brn);
  }

  if (comm_side[ZHI] >= 0)
    MPI_Wait(&req_send,MPI_STATUS_IGNORE);
}

/* ----------------------------------------------------------------------
 * Packing the buffer for reverse communication
 * ----------------------------------------------------------------------*/
void Comm::reverse_pack(Brain *brn, int flag) {
  int i, ag_id;

  int nall = brn->nall;

  double *xhi = brn->xhi;
  double **x = brn->x;

  int *is_loc = brn->is_loc;

  int c = 0;
  send_buf[c++] = ubuf(buf_size[flag]).d;

  for (i=0; i<nall; i++) {
    if (is_loc[i]) continue;
    if (flag == XHI && x[i][0] < xhi[0]) continue;
    else if (flag == YHI && x[i][1] < xhi[1]) continue;
    else if (flag == ZHI && x[i][2] < xhi[2]) continue;

    //printf("proc %i: pack itag = %li, c=%i, i=%i, buf_size =%i \n",
      //     brn->me,brn->tag[i],c,i,buf_size[flag]);

    send_buf[c++] = ubuf(brn->tag[i]).d;
    send_buf[c++] = ubuf(brn->type[i]).d;
    for (ag_id=0; ag_id<num_agents; ag_id++)
      send_buf[c++] = brn->deriv[ag_id][i];
  }

}

/* ----------------------------------------------------------------------
 * Unpacking the buffer for reverse communication
 * ----------------------------------------------------------------------*/
void Comm::reverse_unpack(Brain *brn) {
  int i, ag_id;
  tagint itag;

  int c = 0;
  int unpack_size = (int) ubuf(recv_buf[c++]).i;

  while (c <= unpack_size) {
    itag = (tagint) ubuf(recv_buf[c++]).i;
    i = brn->map[itag];

    //printf("proc %i: unpack itag = %li, c=%i, i=%i, buf_size =%i \n",
         //  brn->me,itag,c-1,i,unpack_size);

    if (i != -1) {
      brn->type[i] = (int) ubuf(recv_buf[c++]).i;
      for (ag_id=0; ag_id<num_agents; ag_id++)
        brn->deriv[ag_id][i] += recv_buf[c++];
    } else
      c += num_agents + 1;

  }

}
