#include "comm_mpi.h"

void border_exchange() {
  MPI_Sendrecv (&lattice[LI(local_grid_width_halo-2, 1)], 1, border_col, east, 0,
                &lattice[LI(0, 1)], 1, border_col, west, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE );

  MPI_Sendrecv (&lattice[LI(1, 1)], 1, border_col, west, 1,
                &lattice[LI(local_grid_width_halo-1, 1)], 1, border_col, east, 1,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE );


  MPI_Sendrecv ( &lattice[LI(0, 1)], 1, border_row, north, 0,
                &lattice[LI(0, local_grid_height_halo-1)], 1, border_row, south, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE );
  MPI_Sendrecv (&lattice[LI(0, local_grid_height_halo-2)], 1, border_row, south, 0,
                &lattice[LI(0, 0)], 1, border_row, north, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
}

void nonb_exchange() {
  
  // EAST
  MPI_Request requestEast, requestWest;

  MPI_Isend(&lattice[LI(local_grid_width_halo-2, 1)]/* east */, 1, border_col, east,
	    0, MPI_COMM_WORLD, &requestEast);
  MPI_Isend(&lattice[LI(1, 1)] /* west */, 1, border_col, west, 
	    1, MPI_COMM_WORLD, &requestWest);
  
  MPI_Irecv(&lattice[LI(0, 1)], 1 /* east */ , border_col, west,
	    0, MPI_COMM_WORLD, &requestWest);
  MPI_Irecv(&lattice[LI(local_grid_width_halo-1, 1)] /* west */, 1, border_col, east,
	    1, MPI_COMM_WORLD, &requestEast);
  
  // Because we are sending corners, sync is needed here
  MPI_Wait(&requestWest, MPI_STATUS_IGNORE);
  MPI_Wait(&requestEast, MPI_STATUS_IGNORE);
  MPI_Barrier(MPI_COMM_WORLD);
  
    
  // NORTH - SOUTH
  MPI_Request requestNorth,  requestSouth;
  MPI_Isend(&lattice[LI(0, 1)] /* north */, 1, border_row, north, 
	    0, MPI_COMM_WORLD, &requestNorth);
  MPI_Isend(&lattice[LI(0, local_grid_height_halo-2)] /* south */ , 1, border_row, south,
	    1, MPI_COMM_WORLD, &requestSouth);
  
  
  MPI_Irecv(&lattice[LI(0, local_grid_height_halo-1)] /* north */, 1, border_row, south,
	    0,  MPI_COMM_WORLD, &requestSouth);
  MPI_Irecv(&lattice[LI(0, 0)] /* south */, 1, border_row, north, 
	    1, MPI_COMM_WORLD,  &requestNorth);
  
  
  MPI_Wait(&requestNorth, MPI_STATUS_IGNORE);
  MPI_Wait(&requestSouth, MPI_STATUS_IGNORE);

  MPI_Barrier(MPI_COMM_WORLD);
  
}

void shmm_exchange() {

  
  MPI_Barrier(nodecomm);
  
   for(int j=0;j < local_grid_height; j++){
     memcpy (&lattice[LI(0, 1+j)], 
	     &westptr[LI(local_grid_width_halo-2, 1+j)], 
	     sizeof(point_t) );
     memcpy (&lattice[LI(local_grid_width_halo-1, 1+j)], 
	     &eastptr[LI(1, 1+j)], 
	     sizeof(point_t) );
   }
  
  MPI_Barrier(nodecomm);
  
  
  for(int i=0; i < local_grid_width_halo; ++i) {
    memcpy ( &lattice[LI(0+i, 0)], 
	     &southptr[LI(0+i, local_grid_height_halo-2)],
	     sizeof(point_t) );
    memcpy ( &lattice[LI(0+i, local_grid_height_halo-1)],
	     &northptr[LI(0+i, 1)], 
	     sizeof(point_t) );
  } 
  // most likely not needed
  MPI_Barrier(nodecomm);
  
}



/* count number of intra and inter node partners */
void get_n_partners ()
{ 
  int j;
  for (j=0; j < n_partners; j++)
    /* If partner has a valid mapping in shm communicator then it is on the same node */
     partners[j] == MPI_UNDEFINED ? n_internode_partners++ : n_intranode_partners++; 
}

/* defines global rank  -> shmcomm rank mapping;
   output: partners_map is array of ranks in shmcomm  */
void translate_ranks()
{
  MPI_Group world_group, shared_group;
  
  /* create MPI groups for global communicator and shm communicator */
  MPI_Comm_group (MPI_COMM_WORLD, &world_group); 
  MPI_Comm_group (nodecomm, &shared_group);
  
  MPI_Group_translate_ranks (world_group, n_partners, partners, shared_group,
			     partners_map); 
}


void get_neighbor_ptrs(){

 MPI_Aint sz;
  int dsp_unit;
  //Get pointer to the north neighbor
  if ( partners_map[0] != MPI_UNDEFINED ){
    MPI_Win_shared_query( nodewin,  partners_map[0], &sz, &dsp_unit,
			  &northptr);
  }
  //Get pointer to the south neighbor
  if ( partners_map[1] != MPI_UNDEFINED) {
    MPI_Win_shared_query( nodewin,  partners_map[1], &sz, &dsp_unit,
			  &southptr);
  }
  //Get pointer to the east neighbor
  if ( partners_map[2] != MPI_UNDEFINED) {
    MPI_Win_shared_query( nodewin,  partners_map[2], &sz, &dsp_unit,
			  &eastptr);
  }
  //Get pointer to the west neighbor
  if ( partners_map[3] != MPI_UNDEFINED){
    MPI_Win_shared_query( nodewin,  partners_map[3], &sz, &dsp_unit,
			  &westptr);
  }
}


/* print number of intra and inter node partners */
void print_n_partners ()
{
  int j, partner;
  char tmp_str_intra[n_partners*16]; 
  char tmp_str_inter[n_partners*10];
  int pos_intra = 0, pos_inter = 0;
  
  for (j=0; j < n_partners; j++)
    {
      partner = partners[j]; /* partner is in the world notation */
      if (partners_map[j] != MPI_UNDEFINED) /* partner j is on the same node  */
	pos_intra += sprintf (&tmp_str_intra[pos_intra], ", %d (%d)", partner,
			      partners_map[j]);
      else
	pos_inter += sprintf (&tmp_str_inter[pos_inter], ", %d", partner);
    }
  
  if (n_internode_partners)
    printf ("I'm rank %d with %d internode partner%c%s \n", 
	    rank, n_internode_partners, n_internode_partners >1?'s':' ', tmp_str_inter);
  
  if (n_intranode_partners) 
    printf ("I'm rank %d with %d intranode partner%c%s\n", 
	    rank, n_intranode_partners, n_intranode_partners > 1?'s':' ', tmp_str_intra);
}


void init_shmmMPI(){

 // Create node-local communicator
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank,
              MPI_INFO_NULL, &nodecomm);
  
  // Check it all went as expected
  MPI_Comm_size(nodecomm, &nodesize);
  MPI_Comm_rank(nodecomm, &noderank);

  printf("Rank %d in COMM_WORLD is rank %d in nodecomm \n",
	 rank, noderank);
  // n_partners allocation space
  partners_map = (int*)malloc(sizeof(int)*n_partners);
  partners     = (int*)malloc(sizeof(int)*n_partners);
  
 /* map neighbors */
  partners[0] = north;
  partners[1] = south;
  partners[2] = east;
  partners[3] = west;
  
  translate_ranks();
  get_n_partners ();
  print_n_partners ();
  // Allocated oldroad as a shared array, contiguous across processes
  winsize = (local_grid_width_halo * local_grid_height_halo )*sizeof(point_t);
  
  // displacements counted in units of integers
  disp_unit = sizeof(point_t);
  MPI_Win_allocate_shared(winsize, disp_unit,
			  MPI_INFO_NULL, nodecomm, &lattice, &nodewin);

  
  get_neighbor_ptrs();

}
