#include "lbm.h"
#include "functions.h"

void
init_directions ( void )
{
  int i;
  for( i=0; i<6; i++ )
    {
      c[i][0] = sin( M_PI * i / 3.0 );
      c[i][1] = cos( M_PI * i / 3.0 );
    }
}


void
init_densities ( void )
{
  int global_x, global_y;
  for ( int y=0; y<local_grid_height; y++ ) {
    for ( int x=0; x<local_grid_width; x++ ) {

      // Store physical coordinates with each point
      global_x = coords[0]*(L_WIDTH/dims[0]) + x;
      global_y = coords[1]*(L_HEIGHT/dims[1]) + y;
      lattice[LIB(x,y)].position[0] =
        (float)(global_x + (global_y&1)*0.5);
      lattice[LIB(x,y)].position[1] = (float) global_y;
      /*
        printf ( "%f %f\n",
        lattice[LIB(x,y)].position[0],
        lattice[LIB(x,y)].position[1]
        );
      */
      for ( int i=0; i<6; i++ ) {
        // Start from even distribution everywhere in non-ghost cells
        if ( ! ghost[GY(y)][GX(x)] ) {
          lattice[LIB(x, y)].density[i][NOW] =
	    lattice[LIB(x, y)].density[i][NEXT] = 1.0/6.0;
        }

        else {
          /* This is a ghost cell */
          int n_x = neighbor_x(GY(y), GX(x), i);
          int n_y = neighbor_y(GY(y), i);

#ifdef DEBUG
          if (n_x >= L_WIDTH+2 || n_y >= L_HEIGHT+2 ||  n_x < 0 || n_y < 0) {
            printf("Illegal n_x,n_y = %d,%d\n", n_x, n_y);
            exit(1);
          }
#endif // DEBUG

          /* Is my neighbor a ghost also? */
          if ( ! ghost[n_y][n_x] ) {
            lattice[LIB(x, y)].density[(i+3)%6][NOW] =
	      lattice[LIB(x, y)].density[i][NEXT] =
	      1.0/6.0; // No
          }

          else {
            lattice[LIB(x, y)].density[(i+3)%6][NOW] =
	      lattice[LIB(x, y)].density[i][NEXT] = 
	      0.0; // Yes
          }
        }
      }
    }
  }
}



void init_mpi() {
  
  /* Get the local rank and the total number of processes */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

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

}


void init_local_sizes() {
  
  local_grid_width = (L_WIDTH / dims[0]);
  local_grid_height = (L_HEIGHT / dims[1]);

  /* Add leftover columns to last rank in row */
  if(coords[0] + 1 == dims[0]) {
    local_grid_width += L_WIDTH % dims[0];
  }
  /* Add leftover rows to last rank in column */
  if(coords[1] + 1== dims[1]) {
    local_grid_height += L_HEIGHT % dims[1];
  }

  /* Local grid sizes with halo  */
  local_grid_width_halo = local_grid_width + 2 * BORDER;
  local_grid_height_halo = local_grid_height + 2 * BORDER;
}

void init_datatypes() {
  
  /* MPI Datatypes */
  MPI_Datatype tmp_type;
  MPI_Type_create_struct(3,
                         (int[]){2,2,2*6},
                         (MPI_Aint[]){offsetof(point_t,position), offsetof(point_t, velocity), 
			     offsetof(point_t, density)},
                         (MPI_Datatype[]){MPI_FLOAT, MPI_FLOAT, MPI_FLOAT},
                         &tmp_type);
  MPI_Aint lb, extent;
  MPI_Type_get_extent(tmp_type, &lb, &extent);
  MPI_Type_create_resized(tmp_type, lb, extent, &mpi_point_t);
  MPI_Type_commit(&mpi_point_t);

  /* Border types */
  MPI_Datatype border_col_tmp;
  MPI_Type_vector(local_grid_height, 1, local_grid_width_halo, mpi_point_t,
		  &border_col_tmp);
  MPI_Type_get_extent(border_col_tmp, &lb, &extent);
  MPI_Type_create_resized(border_col_tmp, lb, extent, &border_col);
  MPI_Type_commit(&border_col);

  MPI_Datatype border_row_tmp;
  MPI_Type_contiguous(local_grid_width_halo, mpi_point_t, &border_row_tmp);
  MPI_Type_get_extent(border_row_tmp, &lb, &extent);
  MPI_Type_create_resized(border_row_tmp, lb, extent, &border_row);
  MPI_Type_commit(&border_row);

  /* Subdomain from local perspective */
  MPI_Type_vector(local_grid_height, local_grid_width, local_grid_width_halo,
                  mpi_point_t, &local_subdomain_t);
  MPI_Type_commit(&local_subdomain_t);
}

void init_datastructures() {
  
  /* TODO: Experiment with calloc here */
  lattice = (point_t*)malloc(local_grid_width_halo * local_grid_height_halo *
			     sizeof(point_t));
  memset(lattice, 0, local_grid_width_halo * local_grid_height_halo *
	 sizeof(point_t));

#ifdef DEBUG
  for (int y = 0; y < local_grid_height_halo; ++y) {
    for (int x = 0; x < local_grid_width_halo; ++x) {
      if (lattice[LI(x, y)].velocity[0] != 0 ||
          lattice[LI(x, y)].velocity[1] != 0 ||
          lattice[LI(x, y)].density[0][NOW] != 0 ||
          lattice[LI(x, y)].density[1][NOW] != 0 ||
          lattice[LI(x, y)].density[2][NOW] != 0 ||
          lattice[LI(x, y)].density[3][NOW] != 0 ||
          lattice[LI(x, y)].density[4][NOW] != 0 ||
          lattice[LI(x, y)].density[5][NOW] != 0 ||
          lattice[LI(x, y)].density[0][NEXT] != 0 ||
          lattice[LI(x, y)].density[1][NEXT] != 0 ||
          lattice[LI(x, y)].density[2][NEXT] != 0 ||
          lattice[LI(x, y)].density[3][NEXT] != 0 ||
          lattice[LI(x, y)].density[4][NEXT] != 0 ||
          lattice[LI(x, y)].density[5][NEXT] != 0) {
        printf("Rank %d: Lattice not zeroed\n", rank);
        exit(1);
      }
    }
  }
#endif // DEBUG
}

