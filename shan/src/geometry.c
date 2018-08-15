#include "lbm.h"
#include "definitions.h"
#include "geometry.h"
#include <stdbool.h>
#include <math.h>

#if WITH_GEOMETRY == Poiseuille
void
init_ghost ( void )
{
    // Base and top lines
    for ( int x=0; x<L_WIDTH; x++ )
    {
        ghost[0][x] = true;
        ghost[L_HEIGHT-1][x] = true;
    }
}
#elif WITH_GEOMETRY == Moffatt
void
init_ghost ( void )
{
    // Base line
    for ( int x=0; x<L_WIDTH; x++ )
        ghost[0][x] = true;

    // Top + wedge
    for ( int x=0; x<(L_WIDTH/2); x++ )
    {
        int stop_y = CHANNEL;
        if ( x>L_WIDTH/4 )
            stop_y += 3*(x-L_WIDTH/4);
        for ( int y=L_HEIGHT-1; y>stop_y; y-- )
        {
            ghost[y][x] = true;
            ghost[y][L_WIDTH-x-1] = true;
        }
    }
}
#elif WITH_GEOMETRY == Cylinder
void
init_ghost ( void )
{
    // Base and top lines
    for ( int x=0; x<L_WIDTH; x++ )
    {
        ghost[0][x] = true;
        ghost[L_HEIGHT-1][x] = true;
    }

    const int
        xmid = SCALE*100, ymid = SCALE*100;
    const float
        radius = SCALE*20;
    for ( int y=0; y<L_HEIGHT; y++ )
    {
        for ( int x=0; x<L_WIDTH; x++ )
        {
            if ( sqrt(pow(x-xmid,2.0)+pow(y-ymid,2.0)) < radius )
                ghost[y][x] = true;
        }
    }

/* Moffatt settings
    // Top + wedge
    for ( int x=0; x<(L_WIDTH/2); x++ )
    {
        int stop_y = CHANNEL;
        if ( x>L_WIDTH/4 )
            stop_y += 3*(x-L_WIDTH/4);
        for ( int y=L_HEIGHT-1; y>stop_y; y-- )
        {
            ghost[y][x] = true;
            ghost[y][L_WIDTH-x-1] = true;
        }
    }
*/
}
#endif
