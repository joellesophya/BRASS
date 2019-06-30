#ifndef HASHES_H
#define HASHES_H

#include "khash.h"
KHASH_MAP_INIT_STR(str,int);

struct HASH {
	khash_t(str) *ind2ind;
};

#endif
