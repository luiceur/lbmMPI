#pragma once
/*
 * Problem-specific parameters for domain geometry and force settings
 */

#define Poiseuille 0
#define Moffatt 1
#define Cylinder 2

#ifndef WITH_GEOMETRY
#define WITH_GEOMETRY Cylinder
#endif

#if WITH_GEOMETRY == Poiseuille
    #ifndef SCALE
    #define SCALE 5
    #endif
    #define L_WIDTH (300*SCALE)
    #define L_HEIGHT (200*SCALE)
    #define MAX_ITER 50000           /* Number of iterations */
    #define CHECKPOINT (MAX_ITER-1)
    #define SNAPSHOT 1000

    #define FORCE 5e-3                /* Magnitude of external force */
    /* Location of external force */
    #define FORCE_COND (GX(x)==1)
    #define FORCE_COND_CUDA (GX_device(x)==1)
#elif WITH_GEOMETRY == Moffatt
    #ifndef SCALE
    #define SCALE 5
    #endif
    #define L_WIDTH (300*SCALE)
    #define L_HEIGHT (300*SCALE)
    #define CHANNEL (L_HEIGHT/4)


    #define MAX_ITER 200000         /* Number of iterations */
    #define CHECKPOINT (MAX_ITER-1)
    #define SNAPSHOT 500

    #define FORCE 1e-3             /* Magnitude of external force */
    /* Location of external force */
    #define FORCE_COND (GX(x)==1 && GY(y)>10 && GY(y)<CHANNEL)
    #define FORCE_COND_CUDA (GX_device(x)==1 && GY_device(y) > 1 && GY_device(y)<CHANNEL)
#elif WITH_GEOMETRY == Cylinder
    #ifndef SCALE
    #define SCALE 1
    #endif
    #define L_WIDTH (300*SCALE)
    #define L_HEIGHT (200*SCALE)
    #define MAX_ITER 20000           /* Number of iterations */
   //#define CHECKPOINT (MAX_ITER-1)
    #define CHECKPOINT 1000
    #define SNAPSHOT 10000

    #define FORCE 5e-3                /* Magnitude of external force */
    /* Location of external force */
    #define FORCE_COND (GX(x)==1)
    #define FORCE_COND_CUDA (GX_device(x)==1)
#endif

void init_ghost ( void );
