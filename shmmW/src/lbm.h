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
extern point_t* lattice;               /* Structure to hold the local grid of densities and velocities */
extern point_t* global_lattice;
extern point_t* lattice_device;

extern bool ghost[L_HEIGHT][L_WIDTH];// = {0};  /* Solid points */
extern bool* ghost_device;

extern float c[6][2];                  /* Neighbor direction vectors */
extern float* c_device;

extern float force[2]; // = {0.0, FORCE};   /* Magnitude of external force */


extern int offsets[2][6][2];
/*
          // Even/odd neighbor coordinate shifts
{
    { {0,1}, {1,0}, {1,-1}, {0,-1}, {-1,-1}, {-1,0} },
    { {0,1}, {1,1}, { 1,0}, {0,-1}, {-1, 0}, {-1,1} }
};
*/

extern int64_t iter;               /* Iteration count */

/* MPI specific globals */
extern int rank;
extern int nprocs;
#define NDIMS (2)               /* Number of dimensions in processor topology */
extern int dims[NDIMS];
extern int* dims_device;
extern int coords[NDIMS];              /* Local coordinates in processor topology */
extern int* coords_device;
extern int north, south, west, east;   /* Cardinal directions in cartesian grid */
extern int local_grid_width;
extern int local_grid_height;
extern int local_grid_width_halo;
extern int local_grid_height_halo;
extern MPI_Comm cartesian; // = MPI_COMM_NULL;      /* Cartesian communicator */
extern MPI_Datatype mpi_point_t;
extern MPI_Datatype local_grid_type_no_borders;
extern MPI_Datatype local_grid_type;
extern MPI_Datatype border_col;
extern MPI_Datatype border_row;
extern MPI_Datatype local_subdomain_t;


/**********************/
/* External functions */
/**********************/

extern void init_device(point_t** lattice_device,
                        point_t* lattice,
                        int local_grid_width_halo,
                        int local_grid_height_halo,
                        size_t sizeof_point_t,
                        bool** ghost_device,
                        bool ghost[][L_WIDTH],
                        float c[][2],
                        float** c_device,
                        int* coords,
                        int** coords_device,
                        int* dims,
                        int** dims_device);

extern float collide( point_t* lattice_device_p,
                         bool* ghost_device_p,
                         float* c_device_p,
                         int* coords_device,
                         int* dims_device,
                         int local_grid_width_device,
                         int local_grid_height_device,
                         int local_grid_width_halo_device,
                         int local_grid_height_halo_device );

extern void copyToHost(point_t* lattice,
                       point_t* lattice_device,
                       int local_grid_width_halo,
                       int local_grid_height_halo,
                       size_t sizeof_point_t);

extern void copyToDevice(point_t* lattice,
                         point_t* lattice_device,
                         int local_grid_width_halo,
                         int local_grid_height_halo,
                         size_t sizeof_point_t);

extern void propagate(point_t* lattice_device_p,
                      int* coords_device,
                      int* dims_device,
                      int local_grid_width_halo_device,
                      int local_grid_height_halo_device,
                      int iter,
                      int rank);

extern void free_device(point_t* lattice_device,
                        float* c_device,
                        bool* ghost_device,
                        int* coords_device,
                        int* dims_device);

extern void host_propagate(void);
extern void host_collide(void);
