#pragma once
/**********************/
/* External variables */
/**********************/

extern point_t* lattice;
extern point_t* global_lattice;
extern point_t* lattice_device;

extern bool ghost[L_HEIGHT][L_WIDTH];
extern bool* ghost_device;

extern float c[6][2];
extern float* c_device;

extern float force[2];

extern int offsets[2][6][2];

extern int64_t iter;

extern int rank;
extern int nprocs;
extern int dims[NDIMS];
extern int* dims_device;
extern int coords[NDIMS];
extern int* coords_device;
extern int north, south, west, east;
extern int local_grid_width;
extern int local_grid_height;
extern int local_grid_width_halo;
extern int local_grid_height_halo;
extern MPI_Comm cartesian;
extern MPI_Datatype mpi_point_t;
extern MPI_Datatype local_grid_type_no_borders;
extern MPI_Datatype local_grid_type;
extern MPI_Datatype border_col;
extern MPI_Datatype border_row;

/**********************/
/* External functions */
/**********************/

extern int local_neighbor_y (int, int);
extern int local_neighbor_x (int, int, int);

inline int local_neighbor_x (int y, int x, int i) {
    if( GY(y-1)%2 ) return ( (x + offsets[1][i][1]));
    else      return ( (x + offsets[0][i][1]));
}

inline int local_neighbor_y (int y, int i) {
    if( GY(y-1)%2 ) return ( (y + offsets[1][i][0]));
    else      return ( (y + offsets[0][i][0]));
}
