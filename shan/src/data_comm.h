#ifndef DATA_COMM_H
#define DATA_COMM_H

#include <mpi.h>

#ifdef USE_SHAN
#include <SHAN_comm.h>
#include <SHAN_type.h>
#endif


typedef struct 
{
  int nProc;
  int iProc;

  int ndomains;
  int ncommdomains;
  int nownpoints;
  int naddpoints;

  /* general comm */
  int *addpoint_owner;
  int *addpoint_id;
  int *commpartner;
  int *sendcount;
  int *recvcount;
  int **recvindex;
  int **sendindex;

  /* comm vars mpi */
  int nreq;
  MPI_Request *req;
  MPI_Status *stat;
  double **recvbuf;
  double **sendbuf;

  /* global stage counter */
  volatile counter_t *recv_counter;
  volatile counter_t *send_counter;

#ifdef USE_SHAN
    shan_neighborhood_t neighborhood_id;
#endif

} comm_data ;


#endif
