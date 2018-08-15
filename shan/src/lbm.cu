#include <stdio.h>
#include <cuda.h>

#include "definitions.h"
#include "geometry.h"

#define GPU_CHECK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{

    if (code != cudaSuccess)
    {
        fprintf(stderr,"GPU ASSERT: %s %d %s %d\n", cudaGetErrorString(code), code, file, line);
        if (abort) exit(code);
    }
}

__constant__ __device__
int offsets_device[2][6][2] =          /* Even/odd neighbor coordinate shifts */
    {
        { {0,1}, {1,0}, {1,-1}, {0,-1}, {-1,-1}, {-1,0} },
        { {0,1}, {1,1}, { 1,0}, {0,-1}, {-1, 0}, {-1,1} }
    };

__constant__ __device__
float force_device[2] = {0.0, FORCE};

__device__
int local_neighbor_x (int y, int x, int i, int* coords_device, int* dims_device) {
    if( GY_device(y-1)%2 ) return ( (x + offsets_device[1][i][1]));
    else      return ( (x + offsets_device[0][i][1]));
}

__device__
int local_neighbor_y (int y, int i, int* coords_device, int* dims_device) {
    if( GY_device(y-1)%2 ) return ( (y + offsets_device[1][i][0]));
    else      return ( (y + offsets_device[0][i][0]));
}

__global__ void __collide ( point_t* lattice_device_p,
                            bool* ghost_device_p,
                            float* c_device_p,
                            int* coords_device,
                            int* dims_device,
                            int local_grid_width_device,
                            int local_grid_height_device,
                            int local_grid_width_halo_device,
                            int local_grid_height_halo_device )
{

    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int x = blockDim.x * blockIdx.x + threadIdx.x;

    if (x >= local_grid_width_device || y >= local_grid_height_device) {
        return;
    }


    float rho = 0.0;
    float uc = 0.0;

#ifdef DEBUG
    if(LIB_device(x,y) >= local_grid_width_halo_device * local_grid_height_halo_device) {
        printf("Illegal x,y = %d,%d\n",x,y);
        printf("LIB(x,y) = %d\n", LIB_device(x,y));
        /* exit(1); */
    }
#endif // DEBUG

    /* Zero the velocity before computing it */
    lattice_device_p[LIB_device(x, y)].velocity[0] = lattice_device_p[LIB_device(x, y)].velocity[1] = 0.0;

    /* Compute velocity unless lattice site is a ghost */
    if ( ! ghost_device_p[GY_device(y) * L_WIDTH + GX_device(x)] )
    {
        for ( int i=0; i<6; i++ )
        {
            rho += lattice_device_p[LIB_device(x, y)].density[i][NOW];
            lattice_device_p[LIB_device(x, y)].velocity[0] += c_device_p[2*i] * lattice_device_p[LIB_device(x, y)].density[i][NOW];
            lattice_device_p[LIB_device(x, y)].velocity[1] += c_device_p[1 + 2*i] * lattice_device_p[LIB_device(x, y)].density[i][NOW];
        }
        /* rho*u = sum_i( Ni*ci ), so divide by rho to find u: */
        lattice_device_p[LIB_device(x, y)].velocity[0] /= rho;
        lattice_device_p[LIB_device(x, y)].velocity[1] /= rho;
    }

    for ( int i=0; i<6; i++ )
    {
        float
            qi_uaub,
            N_eq,
            delta_N;

        float c_device_p_0 = c_device_p[2*i];
        float c_device_p_1 = c_device_p[1 + 2*i];

        qi_uaub =
            ( c_device_p_1 * c_device_p_1 - 0.5 ) * lattice_device_p[LIB_device(x, y)].velocity[1] * lattice_device_p[LIB_device(x, y)].velocity[1] +
            ( c_device_p_1 * c_device_p_0       ) * lattice_device_p[LIB_device(x, y)].velocity[1] * lattice_device_p[LIB_device(x, y)].velocity[0] +
            ( c_device_p_0 * c_device_p_1       ) * lattice_device_p[LIB_device(x, y)].velocity[0] * lattice_device_p[LIB_device(x, y)].velocity[1] +
            ( c_device_p_0 * c_device_p_0 - 0.5 ) * lattice_device_p[LIB_device(x, y)].velocity[0] * lattice_device_p[LIB_device(x, y)].velocity[0];
        uc = lattice_device_p[LIB_device(x, y)].velocity[0] * c_device_p_0 + lattice_device_p[LIB_device(x, y)].velocity[1] * c_device_p_1;

        // Equilibrium, difference
        N_eq = ( rho / 6.0 ) * ( 1.0 + 2.0 * uc + 4.0 * qi_uaub );
        delta_N = LAMBDA * ( lattice_device_p[LIB_device(x, y)].density[i][NOW] - N_eq );

        // Apply external force at boundary
//        if ( GX_device(x)==1 )
        if ( FORCE_COND_CUDA )
            delta_N += (1.0/3.0) * (c_device_p_0*force_device[0] + c_device_p_1*force_device[1]);

        // Reflections at ghosts
        if( ! ghost_device_p[GY_device(y) * L_WIDTH + GX_device(x)] )
            lattice_device_p[LIB_device(x, y)].density[i][NEXT] = lattice_device_p[LIB_device(x, y)].density[i][NOW] + delta_N;
        else
            lattice_device_p[LIB_device(x, y)].density[(i+3)%6][NEXT] = lattice_device_p[LIB_device(x, y)].density[i][NOW];
    }



}

extern "C" void collide( point_t* lattice_device_p,
                         bool* ghost_device_p,
                         float* c_device_p,
                         int* coords_device,
                         int* dims_device,
                         int local_grid_width_device,
                         int local_grid_height_device,
                         int local_grid_width_halo_device,
                         int local_grid_height_halo_device )
{
    dim3 block(BLOCKSIZE, BLOCKSIZE);
    dim3 grid((local_grid_width_halo_device / BLOCKSIZE) + 1, (local_grid_height_halo_device / BLOCKSIZE) + 1);

    __collide<<<grid, block>>>(lattice_device_p,
                               ghost_device_p,
                               c_device_p,
                               coords_device,
                               dims_device,
                               local_grid_width_device,
                               local_grid_height_device,
                               local_grid_width_halo_device,
                               local_grid_height_halo_device );

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        printf("Collide Error: %s\n", cudaGetErrorString(err));
    cudaDeviceSynchronize();

}

__global__
void
__propagate ( point_t* lattice_device_p, int* coords_device, int* dims_device, int local_grid_width_halo_device, int local_grid_height_halo_device, int iter, int rank)
{

    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int x = blockDim.x * blockIdx.x + threadIdx.x;

    if (x >= local_grid_width_halo_device || y >= local_grid_height_halo_device) {
        return;
    }

    for( int i=0; i<6; i++ )
    {
#ifdef DEBUG
        if(GX_device(x) == 0 && GY_device(y) == 0 && iter==0) {
            printf("RANK: %d -- ", rank);
            printf("local neighbor (i=%d), x=%d,y=%d\n",
                   i,
                   local_neighbor_x(y, x, i, coords_device, dims_device),
                   local_neighbor_y(y, i, coords_device, dims_device));
        }
#endif // DEBUG
        int n_x = local_neighbor_x(y, x, i, coords_device, dims_device);
        int n_y = local_neighbor_y(y, i, coords_device, dims_device);

#ifdef DEBUG
        if(n_x < -1 || n_y < -1) {
            printf("Illegal value n_x=%d,n_y=%d\n", n_x, n_y);
            /* exit(1); */
        }
#endif // DEBUG
        if (n_x <= -1 || n_y <= -1 || n_x >= local_grid_width_halo_device || n_y >= local_grid_height_halo_device) {
            continue;
        }

        /* if (ghost[GY(y)][GX(x)]) { */
        /*     continue; */
        /* } */

        /* if (ghost[GY(n_y)][GX(n_x)]) { */
        /*     continue; */
        /* } */

        lattice_device_p[LI_device(n_x, n_y)].density[i][NOW] = lattice_device_p[LI_device(x, y)].density[i][NEXT];



        // Used during debugging: no density should vanish from the system
        //    check_mass();

    }
}

extern "C" void propagate(point_t* lattice_device_p,
                          int* coords_device,
                          int* dims_device,
                          int local_grid_width_halo_device,
                          int local_grid_height_halo_device,
                          int iter,
                          int rank) {
    dim3 block(BLOCKSIZE, BLOCKSIZE);
    dim3 grid((local_grid_width_halo_device / BLOCKSIZE) + 1, (local_grid_height_halo_device / BLOCKSIZE) + 1);
    __propagate<<<grid, block>>>(lattice_device_p,
                                 coords_device,
                                 dims_device,
                                 local_grid_width_halo_device,
                                 local_grid_height_halo_device,
                                 iter,
                                 rank);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        printf("Propagate Error: %s\n", cudaGetErrorString(err));
    cudaDeviceSynchronize();
}

extern "C" void init_device(point_t** lattice_device,
                            point_t* lattice,
                            int local_grid_width_halo,
                            int local_grid_height_halo,
                            size_t sizeof_point_t,
                            bool** ghost_device,
                            bool* ghost,
                            float* c,
                            float** c_device,
                            int* coords,
                            int** coords_device,
                            int* dims,
                            int** dims_device) {
    GPU_CHECK(cudaMalloc(lattice_device, local_grid_width_halo * local_grid_height_halo * sizeof_point_t));
    GPU_CHECK(cudaMemcpy(*lattice_device, lattice, local_grid_width_halo * local_grid_height_halo * sizeof_point_t, cudaMemcpyHostToDevice));

    GPU_CHECK(cudaMalloc(ghost_device, L_HEIGHT * L_WIDTH * sizeof(bool)));
    GPU_CHECK(cudaMemcpy(*ghost_device, ghost, L_HEIGHT*L_WIDTH*sizeof(bool), cudaMemcpyHostToDevice));

    GPU_CHECK(cudaMalloc(c_device, 2 * 6 * sizeof(float)));
    GPU_CHECK(cudaMemcpy(*c_device, c, 2*6*sizeof(float), cudaMemcpyHostToDevice));

    GPU_CHECK(cudaMalloc(coords_device, NDIMS * sizeof(int)));
    GPU_CHECK(cudaMemcpy(*coords_device, coords, NDIMS * sizeof(int), cudaMemcpyHostToDevice));

    GPU_CHECK(cudaMalloc(dims_device, NDIMS * sizeof(int)));
    GPU_CHECK(cudaMemcpy(*dims_device, dims, NDIMS * sizeof(int), cudaMemcpyHostToDevice));

}

extern "C" void copyToHost(point_t* lattice, point_t* lattice_device, int local_grid_width_halo, int local_grid_height_halo, int sizeof_point_t) {
    cudaDeviceSynchronize();
    GPU_CHECK(cudaMemcpy(lattice, lattice_device, local_grid_width_halo * local_grid_height_halo * sizeof_point_t, cudaMemcpyDeviceToHost));
    cudaDeviceSynchronize();
}

extern "C" void copyToDevice(point_t* lattice, point_t* lattice_device, int local_grid_width_halo, int local_grid_height_halo, int sizeof_point_t) {
    cudaDeviceSynchronize();
    GPU_CHECK(cudaMemcpy(lattice_device, lattice, local_grid_width_halo * local_grid_height_halo * sizeof_point_t, cudaMemcpyHostToDevice));
    cudaDeviceSynchronize();
}

extern "C" void free_device(point_t* lattice_device,
                            float* c_device,
                            bool* ghost_device,
                            int* coords_device,
                            int* dims_device) {
    cudaFree(lattice_device);
    cudaFree(c_device);
    cudaFree(ghost_device);
    cudaFree(coords_device);
    cudaFree(dims_device);
}
