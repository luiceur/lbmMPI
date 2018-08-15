#include "lbm.h"


/* Lattice data
 * These structures are global and independent of the solver's
 * configuration (geometry and force)
 */
point_t* lattice;               /* Structure to hold the local grid of densities and velocities */
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
int rank;
int nprocs;
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

/* Private data */
static const char *optstring = "o:";
static char output_dir[256];

/* Private functions */
static void options ( int argc, char **argv );

void print_row(int row, int rank, int iterchk);
void print_status_all(void);
void print_status(int out_rank, int out_nprocs, int out_dims[2], int out_coords[2],
 int out_north, int out_south, int out_east, int out_west,int out_local_grid_width, 
int out_local_grid_height, int out_local_grid_width_halo, int out_local_grid_height_halo);

#define MPI_CHECK(ans) { mpiAssert((ans), __FILE__, __LINE__); }

inline void mpiAssert(int code, const char *file, int line)
{
    if (code != MPI_SUCCESS)
    {
        int length = -1;
        char string[4096];
        MPI_Error_string(code, string, &length);
        printf("MPI ASSERT: %s %s %d\n", string, file, line);
        exit(code);
    }
}


int
neighbor_y ( int y, int i )
{
    if( y%2 ) return ( (y + offsets[1][i][0] + L_HEIGHT) % L_HEIGHT);
    else      return ( (y + offsets[0][i][0] + L_HEIGHT) % L_HEIGHT);
}


int
neighbor_x ( int y, int x, int i )
{
    if( y%2 ) return ( (x + offsets[1][i][1] + L_WIDTH) % L_WIDTH);
    else      return ( (x + offsets[0][i][1] + L_WIDTH) % L_WIDTH);
}


/* Crude approach to local I/O, using token ring + POSIX append */
void
write_local_points ( const char *filename )
{
    // Rank 0 starts the file, others append
    FILE *out = fopen ( filename, (rank==0)?"w":"a" );
    for ( int y=0; y<local_grid_height; y++ )
        for ( int x=0; x<local_grid_width; x++ )
        {
            point_t *p = &(lattice[LIB(x,y)]);
            int written = fwrite ( p, sizeof(point_t), 1, out );
            if ( written != 1 )
                fprintf ( stderr, "Warning, write error at iter %ld, fwrite returned %d\n", 
			  iter, written );
        }
    fclose ( out );
}


void
snapshot_velocity ( const char *filename )
{
    // Rank 0 starts the file, others append
    FILE *out = fopen ( filename, (rank==0)?"w":"a" );
    for ( int y=0; y<local_grid_height; y++ )
        for ( int x=0; x<local_grid_width; x++ )
        {
            float values[4];
            point_t *p = &(lattice[LIB(x,y)]);
            values[0] = p->position[0];
            values[1] = p->position[1];
            values[2] = p->velocity[0];
            values[3] = p->velocity[1];
            int written = fwrite ( values, sizeof(float), 4, out );
            if ( written != 4 )
                fprintf ( stderr, "Warning, write error at iter %ld, fwrite returned %d\n", 
			  iter, written );
        }
    fclose ( out );
}


void
write_checkpoint ( const char *filename, bool complete )
{
    int token = rank, discard, prev, next;
    prev = (rank+nprocs-1)%nprocs;
    next = (rank+1)%nprocs;
    switch ( rank )
    {
        case 0:
            if ( complete )
                write_local_points ( filename );
            else
                snapshot_velocity ( filename );
            if ( nprocs > 1 )
            {
                MPI_Ssend ( &token, 1, MPI_INT, next, 0, MPI_COMM_WORLD );
                MPI_Recv ( &discard, 1, MPI_INT, prev, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE
                );
            }
            break;
        default:
            MPI_Recv ( &discard, 1, MPI_INT, prev, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE
            );
            if ( complete )
                write_local_points ( filename );
            else
                snapshot_velocity ( filename );
            MPI_Ssend ( &token, 1, MPI_INT, next, 0, MPI_COMM_WORLD );
            break;
    }
}


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
                    lattice[LIB(x, y)].density[i][NOW] = lattice[LIB(x, y)].density[i][NEXT] = 1.0/6.0;
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
                        lattice[LIB(x, y)].density[(i+3)%6][NOW] = lattice[LIB(x, y)].density[i][NEXT] = 1.0/6.0; // No
                    }

                    else {
                        lattice[LIB(x, y)].density[(i+3)%6][NOW] = lattice[LIB(x, y)].density[i][NEXT] = 0.0; // Yes
                    }
                }
            }
        }
    }
}




void
check_mass ( void )
{
    float sum = 0.0;
    for( int y=0; y<L_HEIGHT; y++ )
        for( int x=0; x<L_WIDTH; x++ )
            for( int i=0; i<6; i++ )
                sum += lattice[LI(x, y)].density[i][NOW];
    printf ( "mass %e\n", sum );
}




void init_mpi() {
    /* Get the local rank and the total number of processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Create a topology of processors in a 2D grid */
    MPI_Dims_create(nprocs, NDIMS, dims);

    /* Create cartesian communicator */
    int periods[2] = {true, false};
    MPI_Cart_create(MPI_COMM_WORLD, NDIMS, dims, periods, 0, &cartesian);

    /* Determine local process coordinates in topology */
    MPI_Cart_coords(cartesian, rank, NDIMS, coords);

    /* Get shifted source and destination ranks, given a shift direction and amount */
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
                           (MPI_Aint[]){offsetof(point_t,position), offsetof(point_t, velocity), offsetof(point_t, density)},
                           (MPI_Datatype[]){MPI_FLOAT, MPI_FLOAT, MPI_FLOAT},
                           &tmp_type);
    MPI_Aint lb, extent;
    MPI_Type_get_extent(tmp_type, &lb, &extent);
    MPI_Type_create_resized(tmp_type, lb, extent, &mpi_point_t);
    MPI_Type_commit(&mpi_point_t);

    /* Border types */
    MPI_Datatype border_col_tmp;
    MPI_Type_vector(local_grid_height, 1, local_grid_width_halo, mpi_point_t, &border_col_tmp);
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
    lattice = (point_t*)malloc(local_grid_width_halo * local_grid_height_halo * sizeof(point_t));
    memset(lattice, 0, local_grid_width_halo * local_grid_height_halo * sizeof(point_t));

/*
    if(rank == 0) {
        global_lattice = (point_t*)malloc(L_WIDTH * L_HEIGHT * sizeof(point_t));
        memset(global_lattice, 0, L_WIDTH * L_HEIGHT * sizeof(point_t));
    }
*/

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


void border_exchange() {
    MPI_Sendrecv ( 
        &lattice[LI(local_grid_width_halo-2, 1)], 1, border_col, east, 0,
        &lattice[LI(0, 1)], 1, border_col, west, 0, 
        MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );

    MPI_Sendrecv (
        &lattice[LI(1, 1)], 1, border_col, west, 1, 
        &lattice[LI(local_grid_width_halo-1, 1)], 1, border_col, east, 1, 
        MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );
/*
    // EAST
    MPI_Request requestEast;
    MPI_CHECK(MPI_Isend(&lattice[LI(local_grid_width_halo-2, 1)], 1, border_col, east, 0, cartesian, &requestEast));
    MPI_CHECK(MPI_Irecv(&lattice[LI(local_grid_width_halo-1, 1)], 1, border_col, east, 1, cartesian, &requestEast));

    // WEST
    MPI_Request requestWest;
    MPI_CHECK(MPI_Isend(&lattice[LI(1, 1)], 1, border_col, west, 1, cartesian, &requestWest));
    MPI_CHECK(MPI_Irecv(&lattice[LI(0, 1)], 1, border_col, west, 0, cartesian, &requestWest));

    // Because we are sending corners, sync is needed here
    MPI_CHECK(MPI_Wait(&requestWest, MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_Wait(&requestEast, MPI_STATUS_IGNORE));
    MPI_Barrier(MPI_COMM_WORLD);
*/

    MPI_Sendrecv (
        &lattice[LI(0, 1)], 1, border_row, north, 0, 
        &lattice[LI(0, local_grid_height_halo-1)], 1, border_row, south, 0, 
        MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );
    MPI_Sendrecv (
        &lattice[LI(0, local_grid_height_halo-2)], 1, border_row, south, 0, 
        &lattice[LI(0, 0)], 1, border_row, north, 0, 
        MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );
/*
    // NORTH
    MPI_Request requestNorth;
    MPI_CHECK(MPI_Isend(&lattice[LI(0, 1)], 1, border_row, north, 0, cartesian, &requestNorth));
    MPI_CHECK(MPI_Irecv(&lattice[LI(0, 0)], 1, border_row, north, 0, cartesian, &requestNorth));

    // SOUTH
    MPI_Request requestSouth;
    MPI_CHECK(MPI_Isend(&lattice[LI(0, local_grid_height_halo-2)], 1, border_row, south, 0, cartesian, &requestSouth));
    MPI_CHECK(MPI_Irecv(&lattice[LI(0, local_grid_height_halo-1)], 1, border_row, south, 0, cartesian, &requestSouth));

    MPI_CHECK(MPI_Wait(&requestNorth, MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_Wait(&requestSouth, MPI_STATUS_IGNORE));
    MPI_CHECK(MPI_Barrier(cartesian));
*/
}

int
main ( int argc, char **argv )
{
    MPI_Init ( &argc, &argv );
    init_mpi();

    options ( argc, argv );

    init_local_sizes();
    init_datatypes();
    init_datastructures();
    init_directions();
    init_ghost();
    init_densities();

#ifdef DEBUG
    print_status_all();
#endif // DEBUG
    border_exchange();

    do {

        if ( iter>0 && ((iter % CHECKPOINT == 0)||(iter % SNAPSHOT == 0)) )
        {

#ifdef DEBUG
            if (rank == 0)
            {
                printf ( "%ld ", iter );
                fflush(stdout);
            }
#endif // DEBUG
            if ( (iter % CHECKPOINT) == 0 )
            {
                char fname[256];
                sprintf ( fname, "%s/%04ld.dat", output_dir, iter / CHECKPOINT );
                write_checkpoint ( fname, true );
            }
            MPI_Barrier ( MPI_COMM_WORLD );
            if ( (iter % SNAPSHOT) == 0 )
            {
                char fname[256];
                sprintf ( fname, "%s/%04ld.vel", output_dir, iter / SNAPSHOT );
                write_checkpoint ( fname, false );
            }
        }


        host_collide();
	
        border_exchange();

        host_propagate();


    } while ( iter++ < MAX_ITER );

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

    //munmap ( lattice, local_grid_width_halo * local_grid_height_halo * sizeof(point_t) );
    free (lattice);
    MPI_Finalize();

    exit ( EXIT_SUCCESS );
}


static void
options ( int argc, char **argv )
{
    int dir_length = 0;
    memset ( output_dir, 0, 256 );
    if ( rank == 0 )
    {
        int o = 0;
        while ( (o = getopt(argc,argv,optstring)) != -1 )
        {
            switch ( o )
            {
                case 'o':
                    dir_length = strlen ( optarg );
                    strncpy ( output_dir, optarg, dir_length );
                    break;
            }
        }
        if ( dir_length == 0 )
        {
            dir_length = 4;
            strncpy ( output_dir, "data", dir_length );
        }
    }
    MPI_Bcast ( &dir_length, 1, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Bcast ( output_dir, dir_length, MPI_CHAR, 0, MPI_COMM_WORLD );
}


void print_row(int row, int rankchk, int iterchk) {
    if(rank == rankchk && iter == iterchk) {
        for (int x = 0; x < local_grid_width_halo; ++x) {
            printf("%f ", lattice[LI(x,row)].velocity[0]);
        }
        printf("======================================================================\n");
    }
}

void print_status_all(void) {
    if(rank == 0) {
        printf("######################################################################\n");
        printf("#                                STATUS                              #\n");
        printf("######################################################################\n");
    }

    int out_rank = rank;
    int out_nprocs = nprocs;
    int out_dims[2] = {dims[0], dims[1]};
    int out_coords[2] = {coords[0], coords[1]};
    int out_north = north;
    int out_south = south;
    int out_east = east;
    int out_west = west;
    int out_local_grid_width = local_grid_width;
    int out_local_grid_height = local_grid_height;
    int out_local_grid_width_halo = local_grid_width_halo;
    int out_local_grid_height_halo = local_grid_height_halo;

    if(rank == 0) {
        print_status(out_rank, out_nprocs, out_dims, out_coords, out_north, out_south, out_east, out_west, out_local_grid_width, out_local_grid_height, out_local_grid_width_halo, out_local_grid_height_halo);
        for(int i = 1; i < nprocs; i++) {
            MPI_Recv(&out_rank, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&out_nprocs, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(out_dims, 2, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(out_coords, 2, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&out_north, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&out_south, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&out_east, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&out_west, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&out_local_grid_width, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&out_local_grid_height, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&out_local_grid_width_halo, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&out_local_grid_height_halo, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            print_status(out_rank, out_nprocs, out_dims, out_coords, out_north, out_south, out_east, out_west, out_local_grid_width, out_local_grid_height, out_local_grid_width_halo, out_local_grid_height_halo);

        }
    } else {
        MPI_Send(&out_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&out_nprocs, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(out_dims, 2, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(out_coords, 2, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&out_north, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&out_south, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&out_east, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&out_west, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&out_local_grid_width, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&out_local_grid_height, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&out_local_grid_width_halo, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&out_local_grid_height_halo, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

void print_status(int out_rank, int out_nprocs,
                  int out_dims[2], int out_coords[2],
                  int out_north, int out_south, int out_east, int out_west,
                  int out_local_grid_width, int out_local_grid_height,
                  int out_local_grid_width_halo, int out_local_grid_height_halo) {
    printf("RANK: %d\n", out_rank);
    printf("NPROCS: %d\n", out_nprocs);
    printf("DIMENSIONS: (%d,%d)\n", out_dims[0], out_dims[1]);
    printf("COORDS: (%d,%d)\n", out_coords[0], out_coords[1]);
    printf("NORTH, SOUTH, EAST, WEST: %d, %d, %d, %d\n", out_north, out_south, out_east, out_west);
    printf("LOCAL GRID WIDTH, LOCAL GRID HEIGHT: %d, %d\n", out_local_grid_width, out_local_grid_height);
    printf("LOCAL GRID WIDTH + HALO, LOCAL GRID HEIGHT + HALO: %d, %d\n", out_local_grid_width_halo, out_local_grid_height_halo);
    printf("----------------------------------------------------------------------\n");
}
