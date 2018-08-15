typedef float real_t;
typedef int64_t int_t;

#define SCALE (1.0)
#define LAMBDA (-1.0)
#define WIDTH  ((int_t)SCALE*150)
#define HEIGHT ((int_t)SCALE*150)

typedef struct {
    real_t velocity[2];
    real_t density[2][6];
} point_t;

/* Even/odd neighbor coordinate shifts */
int_t offsets[2][6][2] =
{
    { {0,1}, {1,0}, {1,-1}, {0,-1}, {-1,-1}, {-1,0} },  // Even rows
    { {0,1}, {1,1}, { 1,0}, {0,-1}, {-1, 0}, {-1,1} }   // Odd  rows
};

/* Neighbor direction vectors */
real_t c[6][2];

/* Types of points in the domain:
 * - points that are part of the fluid flow
 * - points that are solid, but with a fluid neighbor
 * - points that are solid, with solid neighbors only
 */
typedef enum { FLOW, WALL, SOLID } solid_t;

/* Rectangular section of the domain */
typedef struct {
    int_t origin[2], size[2];
} subdomain_t;

/* Functions found in the main file */
void collide ( subdomain_t *area );
void propagate ( subdomain_t *area );
void allocate_domain ( void );
void free_domain ( void );
void init_directions ( void );
void init_geometry ( void );
void dump_geometry_ppm ( void );
void dump_velocity_vectors ( void );
