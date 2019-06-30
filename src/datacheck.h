#ifndef DATACHECK_H
#define DATACHECK_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "brass.h"

void format_checks(char *pedfilename, struct DATA_STRUCT *data_struct);

#ifdef __cplusplus
}
#endif

#endif
