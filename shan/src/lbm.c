#include "lbm.h"
#include "comm_data.h"
#include "functions.h"
#include "io.h"
#include "comm_mpi.h"
#include "init.h"
#include "options.h"
#include "comm_data.h"

/* Lattice data
 * These structures are global and independent of the solver's
 * configuration (geometry and force)
 */
point_t* lattice;     /* Structure to hold the local grid of densities 
			 and velocities */
point_t* global_lattice;
point_t* lattice_device;

bool ghost[L_HEIGHT][L_WIDTH] = {0};  /* Solid points */
bool* ghost_device;

float c[6][2];                  /* Neighbor direction vectors */
float* c_device;

float force[2] = {0.0, FORCE};   /* Magnitude of external force */

int offsets[2][6][2] =          /* Even/odd neighbor coordinate shifts */
  {
    { {0,1}, {1,0}, {1,-1}, {0,-1}, {-1,-1}, {-1,0} },
    { {0,1}, {1,1}, { 1,0}, {0,-1}, {-1, 0}, {-1,1} }
  };

int64_t iter = 0;               /* Iteration count */

/* MPI specific globals */
int rank,rankMPI;
int nProcMPI, iProcMPI;
int nprocs,sizeMPI;

#define NDIMS (2)               /* Number of dimensions in processor topology */
int dims[NDIMS];
int* dims_device;
int coords[NDIMS];              /* Local coordinates in processor topology */
int* coords_device;
int north, south, west, east;   /* Cardinal directions in cartesian grid */
int local_grid_width;
int local_grid_height;
int local_grid_width_halo;
int local_grid_height_halo;
MPI_Comm cartesian = MPI_COMM_NULL;      /* Cartesian communicator */
MPI_Datatype mpi_point_t;
MPI_Datatype local_grid_type_no_borders;
MPI_Datatype local_grid_type;
MPI_Datatype border_col;
MPI_Datatype border_row;
MPI_Datatype local_subdomain_t;

MPI_Comm nodecomm;
int nodesize, noderank, nodestringlen, upoffset, dnoffset;
MPI_Win nodewin;
MPI_Aint winsize;
int disp_unit;
point_t *northptr, *southptr, *eastptr, *westptr;
int n_partners = 4; //2D
int * partners_map;
int * partners; 
int n_intranode_partners;
int n_internode_partners;

/* Private data */
const char *optstring = "o:";
char output_dir[256];

comm_data* cd;


int
main ( int argc, char **argv )
{

  MPI_Init ( &argc, &argv );
  //
  init_mpi();

  //choose input and parameters
  options ( argc, argv );
  
  init_local_sizes();

  // Init MPI
  init_shmmMPI();

  // Create datatypes for plain MPI
  init_datatypes();
  
  init_directions();
  init_ghost();
  init_densities();
  

  // first halo exchange in plain MPI
  border_exchange();
  
  //MPI epoch 
  MPI_Win_lock_all(0, nodewin);
  MPI_Win_sync(nodewin);
  
  do {
    
    //dump snapshot
    writeIO();
    
    //collide comp
    host_collide();
    
    //shared memory halo exchange
    shmm_exchange();
    
    //propagate comp
    host_propagate();
    

  } while ( iter++ < MAX_ITER );

  //MPI epoch 
  MPI_Win_unlock_all(nodewin);
  
  MPI_Barrier(MPI_COMM_WORLD);

  if(rank==0) {
    free(global_lattice);
    puts(""); // Put newline at end of output
  }


  MPI_Type_free(&mpi_point_t);
  MPI_Type_free(&border_col);
  MPI_Type_free(&border_row);
  MPI_Type_free(&local_subdomain_t);

  MPI_Comm_free(&cartesian);

  //deallocate window memory
  MPI_Win_free(&nodewin);
  
  //free (lattice);
  MPI_Finalize();

  exit ( EXIT_SUCCESS );
}
