#pragma once
/*************/
/* Constants */
/*************/

#define LAMBDA (-1.0)
#define BORDER (1)                /* Thickness of halo */

#define NOW (0)
#define NEXT (1)

#define NDIMS (2)

/* CUDA BLOCK SIZE */
#define BLOCKSIZE 32

/**********/
/* Macros */
/**********/

/* global lattice index */
#define GLI(X, Y) ((X) + L_WIDTH * (Y))
/* local lattice index */
#define LI(X, Y) ((X) + local_grid_width_halo * (Y))
#define LI_device(X, Y) ((X) + local_grid_width_halo_device * (Y))
/* local lattice index skipping border */
#define LIB(X, Y) (((X) + BORDER) + local_grid_width_halo * ((Y) + BORDER))
#define LIB_device(X, Y) (((X) + BORDER) + local_grid_width_halo_device * ((Y) + BORDER))
/* global x */
#define GX(X) ((X) + coords[0] * (L_WIDTH/dims[0]))
#define GX_device(X) ((X) + coords_device[0] * (L_WIDTH/dims_device[0]))
/* global y */
#define GY(Y) ((Y) + coords[1] * (L_HEIGHT/dims[1]))
#define GY_device(Y) ((Y) + coords_device[1] * (L_HEIGHT/dims_device[1]))
/* local x skipping border */
#define LX(X) ((X) + BORDER)
/* local y skipping border */
#define LY(Y) ((Y) + BORDER)

/*********/
/* Types */
/*********/

/* Lattice point data structure
 * Note that this structure is mirrored in an MPI structured type
 * which requires corresponding updates if/when this definition changes.
 */
typedef struct {
    /* Position vectors */
    float position[2];

    /* Velocity vectors */
    float velocity[2];

    /* Density vectors */
    float density[6][2];
} point_t;
