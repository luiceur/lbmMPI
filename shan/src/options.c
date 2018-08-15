#include "lbm.h"

void options ( int argc, char **argv )
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

