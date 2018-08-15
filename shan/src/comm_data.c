#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stddef.h>

#include "comm_data.h"
#include "error_handling.h"

#ifdef USE_SHAN
#ifdef USE_GASPI
#include <GASPI.h>
#endif
#endif




static void init_communication_data(int iProc, int nProc)
{
  ASSERT(cd != NULL);

  /* init comm data */
  cd->nProc = nProc;
  cd->iProc = iProc;

  cd->ndomains = 0;
  cd->ncommdomains = 0;
  cd->nownpoints = 0;
  cd->naddpoints = 0;

  cd->commpartner = NULL;
  cd->sendcount = NULL;
  cd->recvcount = NULL;
  cd->addpoint_owner = NULL;
  cd->addpoint_id = NULL;
  cd->recvindex = NULL;
  cd->sendindex = NULL;

  cd->nreq = 0;
  cd->req = NULL;
  cd->stat = NULL;
  cd->recvbuf = NULL;
  cd->sendbuf = NULL;

  cd->send_counter = NULL;
  cd->recv_counter = NULL;

}

void init_communication(int argc, char *argv[])
{

  /* MPI init */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &nprocs);

  /* Create a topology of processors in a 2D grid */
  MPI_Dims_create(nprocs, NDIMS, dims);
  
  /* Create cartesian communicator */
  int periods[2] = {true, true};
  MPI_Cart_create(MPI_COMM_WORLD, NDIMS, dims, periods, 0, &cartesian);
  
  /* Determine local process coordinates in topology */
  MPI_Cart_coords(cartesian, rank, NDIMS, coords);
  
  /* Get shifted source and destination ranks,
     given a shift direction and amount */
  MPI_Cart_shift(cartesian, 1, 1, &north, &south);
  MPI_Cart_shift(cartesian, 0, 1, &west, &east);
  
#ifdef USE_SHAN
#ifdef USE_GASPI
  gaspi_config_t config;
  SUCCESS_OR_DIE (gaspi_config_get(&config));
  config.build_infrastructure = GASPI_TOPOLOGY_NONE;
  SUCCESS_OR_DIE (gaspi_config_set(config));

  gaspi_rank_t iProc, nProc;
  SUCCESS_OR_DIE (gaspi_proc_init (GASPI_BLOCK));
  SUCCESS_OR_DIE (gaspi_proc_rank (&iProc));
  SUCCESS_OR_DIE (gaspi_proc_num (&nProc));
  ASSERT(iProc == rank);
  ASSERT(nProc == nprocs);
#endif
#endif

  init_communication_data(iProcMPI, nProcMPI, cd);

}

