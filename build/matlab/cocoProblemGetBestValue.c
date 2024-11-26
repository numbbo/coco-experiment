#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

#include "mex.h"

#include "coco.h"

double coco_problem_get_best_value(coco_problem_t *problem);

#define UNPACK_PROBLEM(IDX)                                                                                            \
  if (nrhs < IDX) {                                                                                                    \
    mexErrMsgIdAndTxt(mexFunctionName(), "Missing argument " #IDX ", expected problem.");                              \
  }                                                                                                                    \
  coco_problem_t *problem = (coco_problem_t *)(*((size_t *)mxGetData(prhs[IDX])));                                     \
  if (problem == NULL) {                                                                                               \
    mexErrMsgIdAndTxt(mexFunctionName(), "Invalid problem (NULL) given.");                                             \
  }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  UNPACK_PROBLEM(0);

  /* Return best value objective value */
  double res = coco_problem_get_best_value(problem);
  plhs[0] = mxCreateDoubleScalar(res);
}