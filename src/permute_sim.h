#ifndef PERMUTE_H
#define PERMUTE_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "nrutil.h"
#include "hashes.h"
#include "brass.h"

    struct rank_pair {
        size_t ind;
        double val;
    };

    void permute(struct FAMILY *family, struct HASH hash, struct DATA_STRUCT data_struct, char *filename);
    void print_rep(int i_perm, struct FAMILY *family, struct HASH hash, FILE *file_out);

#ifdef __cplusplus
}
#endif

#endif
