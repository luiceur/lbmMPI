#include "lbm.h"

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
          fprintf ( stderr,
		    "Warning, write error at iter %ld, fwrite returned %d\n",
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


void writeIO(){


  if ( iter>0 && ((iter % CHECKPOINT == 0)||(iter % SNAPSHOT == 0)) )
    {
      
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
}
