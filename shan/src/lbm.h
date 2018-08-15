#pragma once
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <getopt.h>
#include <math.h>
#include <mpi.h>
#ifdef _CRAYC
#include <stddef.h>
#endif

#include <sys/mman.h>

#include "definitions.h"
#include "geometry.h"

/****************/
/* Global state */
/****************/
extern point_t* lattice;  /* Structure to hold the local grid of 
			     densities and velocities */
extern point_t* global_lattice;
extern point_t* lattice_device;

extern bool ghost[L_HEIGHT][L_WIDTH];// = {0};  /* Solid points */
extern bool* ghost_device;

extern float c[6][2];                  /* Neighbor direction vectors */
extern float* c_device;

extern float force[2]; // = {0.0, FORCE};   /* Magnitude of external force */


extern int offsets[2][6][2];

extern int64_t iter;               /* Iteration count */

extern char output_dir[256];
extern const char *optstring;
/* MPI specific globals */
extern int rank;
extern int nprocs;
#define NDIMS (2)        /* Number of dimensions in processor topology */
extern int dims[NDIMS];
extern int* dims_device;
extern int coords[NDIMS]; /* Local coordinates in processor topology */
extern int* coords_device;
extern int north, south, west, east; /* Cardinal directions in cartesian grid */
extern int local_grid_width;
extern int local_grid_height;
extern int local_grid_width_halo;
extern int local_grid_height_halo;
extern MPI_Comm cartesian;      /* Cartesian communicator */
extern MPI_Datatype mpi_point_t;
extern MPI_Datatype local_grid_type_no_borders;
extern MPI_Datatype local_grid_type;
extern MPI_Datatype border_col;
extern MPI_Datatype border_row;
extern MPI_Datatype local_subdomain_t;


/*****  MPI Windows  ****/
extern MPI_Comm nodecomm;
extern int nodesize, noderank, nodestringlen, upoffset, dnoffset;
extern MPI_Win nodewin;
extern MPI_Aint winsize;
extern int disp_unit;
extern point_t *northptr, *southptr, *eastptr, *westptr;
extern int n_partners;
extern int * partners_map;
extern int * partners; 
extern int n_intranode_partners;
extern int n_internode_partners;

/**********************/
/* External functions */
/**********************/


extern void host_propagate(void);
extern void host_collide(void);
