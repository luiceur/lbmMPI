#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <mpi.h>

#include "definitions.h"
#include "geometry.h"
#include "host_kernel.h"

void
host_propagate ( void )
{
  for( int y=0; y<local_grid_height_halo; y++ )
    for( int x=0; x<local_grid_width_halo; x++ )
      {
        for( int i=0; i<6; i++ )
          {
#ifdef DEBUG
            if(GX(x) == 0 && GY(y) == 0 && iter==0) {
              printf("RANK: %d -- ", rank);
              printf("local neighbor (i=%d), x=%d,y=%d\n",
                     i,
                     local_neighbor_x(y, x, i),
                     local_neighbor_y(y, i));
            }
#endif // DEBUG
            int n_x = local_neighbor_x(y, x, i);
            int n_y = local_neighbor_y(y, i);

#ifdef DEBUG
            if(n_x < -1 || n_y < -1) {
              printf("Illegal value n_x=%d,n_y=%d\n", n_x, n_y);
              exit(1);
            }
#endif // DEBUG
            if (n_x <= -1 || n_y <= -1 || n_x >= local_grid_width_halo || 
		n_y >= local_grid_height_halo) {
              continue;
            }

            lattice[LI(n_x, n_y)].density[i][NOW] = lattice[LI(x, y)].density[i][NEXT];
          }
      }


}

void
host_collide ( void )
{
#ifdef DATAP
#pragma omp parallel for
#endif // DATAP
  for ( int y=0; y<local_grid_height; y++ )
    {
#ifdef TASKP
#pragma omp task
#endif // TASKP
      for ( int x=0; x<local_grid_width; x++ )
        {
          float rho = 0.0;
          float uc = 0.0;

#ifdef DEBUG
          if(LIB(x,y) >= local_grid_width_halo * local_grid_height_halo) {
            printf("Illegal x,y = %d,%d\n",x,y);
            printf("LIB(x,y) = %d\n", LIB(x,y));
            exit(1);
          }
#endif // DEBUG

          /* Zero the velocity before computing it */
          lattice[LIB(x, y)].velocity[0] = lattice[LIB(x, y)].velocity[1] = 0.0;

          /* Compute velocity unless lattice site is a ghost */
          if ( ! ghost[GY(y)][GX(x)] )
            {
              for ( int i=0; i<6; i++ )
                {
                  rho += lattice[LIB(x, y)].density[i][NOW];
                  lattice[LIB(x, y)].velocity[0] += c[i][0] * lattice[LIB(x, y)].density[i][NOW];
                  lattice[LIB(x, y)].velocity[1] += c[i][1] * lattice[LIB(x, y)].density[i][NOW];
                }
              /* rho*u = sum_i( Ni*ci ), so divide by rho to find u: */
              lattice[LIB(x, y)].velocity[0] /= rho;
              lattice[LIB(x, y)].velocity[1] /= rho;
            }

          for ( int i=0; i<6; i++ )
            {
              float
                qi_uaub,
                N_eq,
                delta_N;
              qi_uaub =
                ( c[i][1] * c[i][1] - 0.5 ) * lattice[LIB(x, y)].velocity[1] * 
		lattice[LIB(x, y)].velocity[1] +
                ( c[i][1] * c[i][0]       ) * lattice[LIB(x, y)].velocity[1] * 
		lattice[LIB(x, y)].velocity[0] +
                ( c[i][0] * c[i][1]       ) * lattice[LIB(x, y)].velocity[0] * 
		lattice[LIB(x, y)].velocity[1] +
                ( c[i][0] * c[i][0] - 0.5 ) * lattice[LIB(x, y)].velocity[0] * 
		lattice[LIB(x, y)].velocity[0];
              uc = lattice[LIB(x, y)].velocity[0] * c[i][0] + lattice[LIB(x, y)].velocity[1] * 
		c[i][1];

              // Equilibrium, difference
              N_eq = ( rho / 6.0 ) * ( 1.0 + 2.0 * uc + 4.0 * qi_uaub );
              delta_N = LAMBDA * ( lattice[LIB(x, y)].density[i][NOW] - N_eq );

              // Apply external force at one end of the domain
              if ( FORCE_COND )
                delta_N += (1.0/3.0) * (c[i][0]*force[0] + c[i][1]*force[1]);

              // Reflections at ghosts
              if( ! ghost[GY(y)][GX(x)] )
                lattice[LIB(x, y)].density[i][NEXT] = lattice[LIB(x, y)].density[i][NOW] + delta_N;
              else
                lattice[LIB(x, y)].density[(i+3)%6][NEXT] = lattice[LIB(x, y)].density[i][NOW];
            }

        }
    }
}
