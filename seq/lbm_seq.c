#define _XOPEN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "lbm_seq.h"

int_t size[2] = { WIDTH, HEIGHT };
point_t *domain;
#define DOM(i,j) domain[(i)*WIDTH+(j)]
#define NOW 0
#define NXT 1

solid_t *geometry;
#define GEOM(i,j) geometry[(i)*WIDTH+(j)]


void
collide ( subdomain_t *area )
{
    int_t
        oy = area->origin[0],
        ox = area->origin[1],
        my = area->origin[0] + area->size[0],
        mx = area->origin[1] + area->size[1];

    for ( int_t y=oy; y<my; y++ )
    {
        for ( int_t x=ox; x<mx; x++ )
        {
            /* Skip computation for solid interior */
            if ( GEOM(y,x) == FLOW || GEOM(y,x) == WALL )
            {
                /* Zero initial values before computation */
                real_t rho = 0.0;
                DOM(y,x).velocity[0] = DOM(y,x).velocity[1] = 0.0;

                /* Compute velocity in flow */
                if ( GEOM(y,x) == FLOW )
                {
                    for ( int_t n=0; n<6; n++ )
                    {
                        rho += DOM(y,x).density[NOW][n];
                        DOM(y,x).velocity[0] +=
                            c[n][0] * DOM(y,x).density[NOW][n];
                        DOM(y,x).velocity[1] +=
                            c[n][1] * DOM(y,x).density[NOW][n];
                    }
                    DOM(y,x).velocity[0] /= rho;
                    DOM(y,x).velocity[1] /= rho;
                }

                /* Compute equilibrium transfer of density */
                for ( int_t n=0; n<6; n++ )
                {
                    real_t qi_uaub, uc, N_eq, delta_N;
                    qi_uaub = 
                        ( c[n][1] * c[n][1] - 0.5 ) * DOM(y,x).velocity[1] * DOM(y,x).velocity[1] +
                        ( c[n][1] * c[n][0]       ) * DOM(y,x).velocity[1] * DOM(y,x).velocity[0] +
                        ( c[n][0] * c[n][1]       ) * DOM(y,x).velocity[0] * DOM(y,x).velocity[1] +
                        ( c[n][0] * c[n][0] - 0.5 ) * DOM(y,x).velocity[0] * DOM(y,x).velocity[0];
                    uc = DOM(y,x).velocity[0] * c[n][0] + DOM(y,x).velocity[1] * c[n][1];
                    N_eq = (rho/6.0) * (1.0 + 2.0 * uc + 4.0 * qi_uaub);
                    delta_N = LAMBDA * (DOM(y,x).density[NOW][n] - N_eq);

                    /* Apply external force at horizontal boundary */
                    if ( x == 0 )
                        delta_N += (1.0/3.0) * ( c[n][0]*0.0 + c[n][1]*5e-3 );

                    if ( GEOM(y,x) == FLOW )        /* Propagate density in flow */
                        DOM(y,x).density[NXT][n] = DOM(y,x).density[NOW][n] + delta_N;
                    else if ( GEOM(y,x) == WALL )   /* Reflect density at walls */
                        DOM(y,x).density[NXT][(n+3)%6] = DOM(y,x).density[NOW][n];
                }
            }
        }
    }
}


void
propagate ( subdomain_t *area )
{
    int_t
        oy = area->origin[0],
        ox = area->origin[1],
        my = area->origin[0] + area->size[0],
        mx = area->origin[1] + area->size[1];
    for ( int_t y=oy; y<my; y++ )
    {
        for ( int_t x=ox; x<mx; x++ )
        {
            for ( int_t n=0; n<6; n++ )
            {
                /* Find neighbor coordinates */
                int_t 
                    ny = y + offsets[(y&1)][n][0],
                    nx = x + offsets[(y&1)][n][1];
                /* Periodic boundaries */
                ny = (ny+HEIGHT)%HEIGHT;
                nx = (nx+WIDTH)%WIDTH;
                /* Update densities to those computed in collision stage */
                DOM(ny,nx).density[NOW][n] = DOM(y,x).density[NXT][n];
            }
        }
    }
}


int
main ( int argc, char **argv )
{
    int_t max_iter = 200;
    if ( argc > 1 )
        max_iter = strtol ( argv[1], NULL, 10 );

    allocate_domain();
    init_geometry();
    init_directions();

    subdomain_t *all_points = malloc ( sizeof(subdomain_t) );
    all_points->origin[0] = all_points->origin[1] = 0;
    all_points->size[0] = HEIGHT;
    all_points->size[1] = WIDTH;

    for ( int_t iter=0; iter<max_iter; iter++ )
    {
        collide ( all_points );
        propagate ( all_points );
    }

    dump_geometry_ppm();
    dump_velocity_vectors();

    free ( all_points );
    free_domain();

    exit ( EXIT_SUCCESS );
}


void
allocate_domain ( void )
{
    domain = (point_t *) malloc ( size[0]*size[1]*sizeof(point_t) );
    geometry = (solid_t *) malloc ( size[0]*size[1]*sizeof(solid_t) );
    memset ( domain, 0, size[0]*size[1]*sizeof(point_t) );
    memset ( geometry, 0, size[0]*size[1]*sizeof(solid_t) );
}


void
free_domain ( void )
{
    free (domain);
    free (geometry);
}


void
init_directions ( void )
{
    for( int_t i=0; i<6; i++ )
    {
        c[i][0] = sin( M_PI * i / 3.0 );
        c[i][1] = cos( M_PI * i / 3.0 );
    }
}


void
init_geometry ( void )
{
    /* Walls at vertical boundaries */
    for ( int_t x=0; x<WIDTH; x++ )
    {
        GEOM(0,x) = SOLID;
        GEOM(HEIGHT-1,x) = SOLID;
    }

    /* Body of solids in the domain */

    /* Sinusoidal wall */
    for ( int_t y=0; y<HEIGHT; y++ )
        for ( int_t x=0; x<WIDTH; x++ )
        {
            if ( y < (HEIGHT/2)+(HEIGHT/3)*sin(2*M_PI*x/(real_t)WIDTH) )
                GEOM(y,x) = SOLID;
        }

    /* Square box in the middle
    for ( int_t y=HEIGHT/4; y<3*HEIGHT/4; y++ )
        for ( int_t x=WIDTH/4; x<3*WIDTH/4; x++ )
            GEOM(y,x) = SOLID;
    */

    /* Detect solids that are walls, i.e. those with neighbors in the flow */
    for ( int_t y=0; y<HEIGHT; y++ )
        for ( int_t x=0; x<WIDTH; x++ )
        {
            if ( GEOM(y,x) == SOLID )
            {
                for ( int_t neigh=0; neigh<6; neigh++ )
                {
                    /* Find neighbor coordinates */
                    int_t 
                        ny = y + offsets[(y&1)][neigh][0],
                        nx = x + offsets[(y&1)][neigh][1];
                    /* Periodic boundaries */
                    ny = (ny+HEIGHT)%HEIGHT;
                    nx = (nx+WIDTH)%WIDTH;
                    /* Set solid to be wall if this neighbor is in the flow */
                    if ( GEOM(ny,nx) == FLOW )
                        GEOM(y,x) = WALL;
                }
            }
        }

    /* Initial density distribution at equilibrium */
    for ( int_t y=0; y<HEIGHT; y++ )
        for ( int_t x=0; x<WIDTH; x++ )
            for ( int_t n=0; n<6; n++ )
                DOM(y,x).density[NOW][n] = 1.0 / 6.0;
}


void
dump_geometry_ppm ( void )
{
    uint8_t triplet[3];
    FILE *out = fopen ( "geometry.ppm", "w" );
    fprintf ( out, "P6 %ld %ld 255\n", WIDTH, HEIGHT );
    for ( int_t y=0; y<HEIGHT; y++ )
    {
        for ( int_t x=0; x<WIDTH; x++ )
        {
            switch ( GEOM(y,x) )
            {
                case FLOW: triplet[0]=triplet[1]=triplet[2] = 0;        break;
                case WALL: triplet[0]=255; triplet[1] = triplet[2] = 0; break;
                case SOLID: triplet[0]=triplet[1]=triplet[2]=255;       break;
            }
            fwrite ( triplet, sizeof(uint8_t), 3, out );
        }
    }
    fclose (out);
}


void
dump_velocity_vectors ( void )
{
    FILE *out = fopen ( "vectors.txt", "w" );
    for ( int_t y=0; y<HEIGHT; y++ )
    {
        for ( int_t x=0; x<WIDTH; x++ )
        {
            real_t norm = sqrt (
                pow(DOM(y,x).velocity[0],2.0) + pow(DOM(y,x).velocity[1],2.0)
            );
            fprintf ( out, "%e %e %e %e\n", x + 0.5*(y&1), (real_t)y,
                DOM(y,x).velocity[1]/norm, DOM(y,x).velocity[0]/norm
            );
        }
    }
    fclose ( out );
}
