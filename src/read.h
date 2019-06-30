#ifndef READ_H
#define READ_H


#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nrutil.h"
#include "hashes.h"
#include "brass.h"

void readpheno(char* filename, struct FAMILY *family, struct HASH *hash, struct DATA_STRUCT data_struct);
void readeigdecomp(char *filename, struct FAMILY *family);
void readkin(char *filename, struct FAMILY *family, struct HASH hash);

#ifdef __cplusplus
}
#endif

#endif
