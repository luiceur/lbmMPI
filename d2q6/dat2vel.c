#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "src/definitions.h"
#include "src/geometry.h"


int
main ( int argc, char **argv )
{
    if ( argc < 2 )
    {
        fputs ( "Need a file name\n", stderr );
        exit ( EXIT_FAILURE );
    }
    point_t *p = malloc ( sizeof(point_t) );

    FILE *in = fopen ( argv[1], "r" );
    for ( int i=0; i<L_WIDTH*L_HEIGHT; i++ )
    {
        int read = fread ( p, sizeof(point_t), 1, in );
        if ( read != 1 )
            fprintf ( stderr, "Warning, read error at point %d\n", i );
        float
            x = p->position[0],
            y = p->position[1],
            vx = p->velocity[0],
            vy = p->velocity[1];
            float norm = sqrt((vx*vx)+(vy*vy));
            if ( norm < 1e-8 )
                norm = 1e-8;
            printf ( "%f %f %e %e\n", x, y, vy/norm, vx/norm );
    }
    fclose ( in );
    free ( p );
    exit ( EXIT_SUCCESS );
}
