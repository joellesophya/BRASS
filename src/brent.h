#ifndef BRENT_H
#define BRENT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#include "compute_residuals.h"

double brent(function_xi *f, double *root);

#ifdef __cplusplus
}
#endif

#endif
