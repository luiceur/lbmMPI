#include "lbm.h"

void border_exchange();
void nonb_exchange();
void shmm_exchange();
void get_n_partners();
void translate_ranks();
void get_neighbor_ptrs();
void print_n_partners();
void init_shmmMPI();
