#include "lbm.h"


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



