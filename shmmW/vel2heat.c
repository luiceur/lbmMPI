#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "src/definitions.h"
#include "src/geometry.h"


float *mat, *mat2;
#define MAT(y,x) mat[(y)*L_WIDTH+(x)]
#define MAT2(y,x) mat2[(y)*L_WIDTH+(x)]


int
main ( int argc, char **argv )
{
    if ( argc < 2 )
    {
        fputs ( "Need a file name\n", stderr );
        exit ( EXIT_FAILURE );
    }

    mat  = malloc ( L_WIDTH*L_HEIGHT*sizeof(float) );

    int length = strlen(argv[1]);
    char fname[length+1];
    memcpy ( fname, argv[1], length );
    fname[length] = 0;

    fname[length-4] = '.';
    fname[length-3] = 'm';
    fname[length-2] = 'a';
    fname[length-1] = 't';

    /* Populate the matrix with z-values
     * Using velocity vector norm to highlight areas with fast flow
     */
    FILE *in = fopen ( argv[1], "r" );
    for ( int i=0; i<(L_HEIGHT*L_WIDTH); i++ )
    {
        point_t *p, point;
        p = &point;
        int read;
        read = fread ( &(p->position[0]), sizeof(float), 1, in );
        read = fread ( &(p->position[1]), sizeof(float), 1, in );
        read = fread ( &(p->velocity[0]), sizeof(float), 1, in );
        read = fread ( &(p->velocity[1]), sizeof(float), 1, in );
        if ( read != 1 )
            fprintf ( stderr, "Warning, read error at %d\n", i );
        float z = 0.0;
        z = sqrtf (powf(p->velocity[0],2.0) + powf(p->velocity[1],2.0));
        int x = (int) floor(p->position[0]);
        int y = (int) floor(p->position[1]);
        MAT(y,x) = z;
    }
    fclose ( in );

    /* Write to gnuplot binary matrix file */
    FILE *out = fopen ( fname, "w" );

    float tmp = (float)L_WIDTH;
    fwrite ( &tmp, sizeof(float), 1, out );
    for ( int j=0; j<L_WIDTH; j++ )
    {
        tmp = (float)j;
        fwrite ( &tmp, sizeof(float), 1, out );
    }
    for ( int i=0; i<L_HEIGHT; i++ )
    {
        tmp = (float)i;
        fwrite ( &tmp, sizeof(float), 1, out );
        for ( int j=0; j<L_WIDTH; j++ )
            fwrite ( &MAT(i,j), sizeof(float), 1, out );
    }

    fclose ( out );
    free ( mat );
    free ( mat2 );
    exit ( EXIT_SUCCESS );
}
