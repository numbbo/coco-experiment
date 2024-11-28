#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <stdlib.h>

#include "coco.cpp"

#include "mex.h"

#define UNPACK_THING(TYPE, THING, IDX)                                                                                 \
  if (nrhs < IDX) {                                                                                                    \
    mexErrMsgIdAndTxt(mexFunctionName(), "Missing argument " #IDX ", expected " #THING ".");                           \
  }                                                                                                                    \
  if (!mxIsScalar(prhs[IDX])) {                                                                                        \
    mexErrMsgIdAndTxt(mexFunctionName(), "Argument " #IDX " not scalar. Expected pointer to " #THING ".");             \
  }                                                                                                                    \
  TYPE *THING = reinterpret_cast<TYPE *>((uint64_t)mxGetScalar(prhs[IDX]));

#define UNPACK_THING_NOT_NULL(TYPE, THING, IDX)                                                                        \
  UNPACK_THING(TYPE, THING, IDX)                                                                                       \
  if (THING == nullptr) {                                                                                                 \
    mexErrMsgIdAndTxt(mexFunctionName(), "Invalid " #THING " (NULL) passed in.");                                      \
  }

#define UNPACK_PROBLEM(IDX) UNPACK_THING_NOT_NULL(coco_problem_t, problem, IDX)
#define UNPACK_OBSERVER(IDX) UNPACK_THING_NOT_NULL(coco_observer_t, observer, IDX)
#define UNPACK_SUITE(IDX) UNPACK_THING_NOT_NULL(coco_suite_t, suite, IDX)

#define RETURN_POINTER(PTR, IDX) plhs[IDX] = mxCreateDoubleScalar(reinterpret_cast<int64_t>(PTR));

void cocoEvaluateFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* const char *class_name = NULL; */
  int nb_objectives;

  /* check for proper number of arguments */
  if (nrhs != 2) {
    mexErrMsgIdAndTxt("cocoEvaluateFunction:nrhs", "Two inputs required.");
  }
  UNPACK_PROBLEM(0);

  /* make sure the second input argument is array of doubles */
  if (!mxIsDouble(prhs[1])) {
    mexErrMsgIdAndTxt("cocoEvaluateFunction:notDoubleArray", "Input x must be an array of doubles.");
  }
  /* test if input dimension is consistent with problem dimension */
  if (!(mxGetN(prhs[1]) == 1 && mxGetM(prhs[1]) == coco_problem_get_dimension(problem)) &&
      !(mxGetM(prhs[1]) == 1 && mxGetN(prhs[1]) == coco_problem_get_dimension(problem))) {
    mexErrMsgIdAndTxt("cocoEvaluateFunction:wrongDimension", "Input x does not comply with problem dimension.");
  }
  /* get the x vector */
  double *x = mxGetPr(prhs[1]);
  /* prepare the return value */
  nb_objectives = coco_problem_get_number_of_objectives(problem);
  plhs[0] = mxCreateDoubleMatrix(1, (size_t)nb_objectives, mxREAL);
  double *y = mxGetPr(plhs[0]);
  /* call coco_evaluate_function(...) */
  coco_evaluate_function(problem, x, y);
}

void cocoEvaluateConstraint(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* const char *class_name = NULL; */
  int nb_constraints;
  double *x;
  double *y;

  /* check for proper number of arguments */
  if (nrhs != 2) {
    mexErrMsgIdAndTxt("cocoEvaluateConstraint:nrhs", "Two inputs required.");
  }
  /* get the problem */
  UNPACK_PROBLEM(0);

  /* make sure the second input argument is array of doubles */
  if (!mxIsDouble(prhs[1])) {
    mexErrMsgIdAndTxt("cocoEvaluateConstraint:notDoubleArray", "Input x must be an array of doubles.");
  }
  /* test if input dimension is consistent with problem dimension */
  if (!(mxGetN(prhs[1]) == 1 && mxGetM(prhs[1]) == coco_problem_get_dimension(problem)) &&
      !(mxGetM(prhs[1]) == 1 && mxGetN(prhs[1]) == coco_problem_get_dimension(problem))) {
    mexErrMsgIdAndTxt("cocoEvaluateConstraint:wrongDimension", "Input x does not comply with problem dimension.");
  }
  /* get the x vector */
  x = mxGetPr(prhs[1]);
  /* prepare the return value */
  nb_constraints = coco_problem_get_number_of_constraints(problem);
  plhs[0] = mxCreateDoubleMatrix(1, (size_t)nb_constraints, mxREAL);
  y = mxGetPr(plhs[0]);
  /* call coco_evaluate_constraint(...) */
  coco_evaluate_constraint(problem, x, y);
}

void cocoRecommendSolution(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check for proper number of arguments */
  if (nrhs != 2) {
    mexErrMsgIdAndTxt("cocoRecommendSolution:nrhs", "Two inputs required.");
  }
  /* get the problem */
  UNPACK_PROBLEM(0);

  /* make sure the second input argument is array of doubles */
  if (!mxIsDouble(prhs[1])) {
    mexErrMsgIdAndTxt("cocoRecommendSolution:notDoubleArray", "Input x must be an array of doubles.");
  }
  /* test if input dimension is consistent with problem dimension */
  if (!(mxGetN(prhs[1]) == 1 && mxGetM(prhs[1]) == coco_problem_get_dimension(problem)) &&
      !(mxGetM(prhs[1]) == 1 && mxGetN(prhs[1]) == coco_problem_get_dimension(problem))) {
    mexErrMsgIdAndTxt("cocoRecommendSolution:wrongDimension", "Input x does not comply with problem dimension.");
  }
  /* get the x vector */
  double *x = mxGetPr(prhs[1]);
  /* call coco_recommend_solution(...) */
  coco_recommend_solution(problem, x);
}

void cocoObserver(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  coco_observer_t *observer = nullptr;

  /* check for proper number of arguments */
  if (nrhs != 2) {
    mexErrMsgIdAndTxt("cocoObserver:nrhs", "Two inputs required.");
  }
  /* get observer_name */
  char *observer_name = mxArrayToString(prhs[0]);
  /* get observer_options */
  char *observer_options = mxArrayToString(prhs[1]);
  /* call coco_observer() */
  observer = coco_observer(observer_name, observer_options);
  /* prepare the return value */
  RETURN_POINTER(observer, 0);
}

void cocoObserverFree(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  UNPACK_OBSERVER(0);

  /* call coco_observer_free() */
  coco_observer_free(observer);
}

void cocoObserverSignalRestart(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* get the observer */
  UNPACK_OBSERVER(0);
  UNPACK_PROBLEM(1);

  /* call coco_observer_signal_restart() */
  coco_observer_signal_restart(observer, problem);
}

void cocoProblemAddObserver(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  coco_problem_t *observedproblem;

  /* check for proper number of arguments */
  if (nrhs != 2) {
    mexErrMsgIdAndTxt("cocoProblemAddObserver:nrhs", "Two inputs required.");
  }
  
  UNPACK_PROBLEM(0);
  UNPACK_OBSERVER(1);

  /* call coco_problem_add_observer() */
  observedproblem = coco_problem_add_observer(problem, observer);

  /* prepare the return value */
  RETURN_POINTER(observedproblem, 0);
}

void cocoProblemFinalTargetHit(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  const mwSize dims[2] = {1, 1};
  size_t *res;

  /* check for proper number of arguments */
  if (nrhs != 1) {
    mexErrMsgIdAndTxt("cocoProblemFinalTargetHit:nrhs", "One input required.");
  }
  /* get the problem */
  UNPACK_PROBLEM(0);

  /* prepare the return value */
  plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
  res = (size_t *)mxGetData(plhs[0]);
  res[0] = coco_problem_final_target_hit(problem);
}

void cocoProblemFree(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  UNPACK_PROBLEM(0);
  /* call coco_problem_free() */
  coco_problem_free(problem);
}

void cocoProblemGetDimension(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  const mwSize dims[2] = {1, 1};
  size_t *res;

  UNPACK_PROBLEM(0);

  /* prepare the return value */
  plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
  res = (size_t *)mxGetData(plhs[0]);
  res[0] = coco_problem_get_dimension(problem);
}

void cocoProblemGetEvaluations(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  const mwSize dims[2] = {1, 1};
  size_t *res;

  /* check for proper number of arguments */
  if (nrhs != 1) {
    mexErrMsgIdAndTxt("cocoProblemGetEvaluations:nrhs", "One input required.");
  }
  UNPACK_PROBLEM(0);

  /* prepare the return value */
  plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
  res = (size_t *)mxGetData(plhs[0]);
  res[0] = coco_problem_get_evaluations(problem);
}

void cocoProblemGetEvaluationsConstraints(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  const mwSize dims[2] = {1, 1};
  size_t *res;

  UNPACK_PROBLEM(0);

  /* prepare the return value */
  plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
  res = (size_t *)mxGetData(plhs[0]);
  res[0] = coco_problem_get_evaluations_constraints(problem);
}

void cocoProblemGetId(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  const char *res;

  UNPACK_PROBLEM(0);
  /* call coco_problem_get_id(...) */
  res = coco_problem_get_id(problem);
  /* prepare the return value */
  plhs[0] = mxCreateString(res);
}

void cocoProblemGetInitialSolution(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  size_t nb_dim;
  double *v; /* intermediate variable that allows to set plhs[0] */

  UNPACK_PROBLEM(0);
  nb_dim = coco_problem_get_dimension(problem);
  plhs[0] = mxCreateDoubleMatrix(1, nb_dim, mxREAL);
  v = mxGetPr(plhs[0]);
  /* call coco_problem_get_initial_solution(...) */
  coco_problem_get_initial_solution(problem, v);
}

void cocoProblemGetLargestFValuesOfInterest(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  UNPACK_PROBLEM(0);

  size_t nb_obj = coco_problem_get_number_of_objectives(problem);
  plhs[0] = mxCreateDoubleMatrix(1, nb_obj, mxREAL);
  double *v = mxGetPr(plhs[0]);
  /* call coco_problem_get_largest_fvalues_of_interest(...) */
  double const *res = coco_problem_get_largest_fvalues_of_interest(problem);
  for (size_t i = 0; i < nb_obj; i++) {
    v[i] = res[i];
  }
}

void cocoProblemGetLargestValuesOfInterest(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  size_t nb_dim;

  UNPACK_PROBLEM(0);

  nb_dim = coco_problem_get_dimension(problem);
  plhs[0] = mxCreateDoubleMatrix(1, nb_dim, mxREAL);
  double *v = mxGetPr(plhs[0]);
  /* call coco_problem_get_largest_values_of_interest(...) */
  double const *res = coco_problem_get_largest_values_of_interest(problem);
  for (size_t i = 0; i < nb_dim; i++) {
    v[i] = res[i];
  }
}

void cocoProblemGetName(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  size_t *ref;
  coco_problem_t *pb = nullptr;
  const char *res;

  /* check for proper number of arguments */
  if (nrhs != 1) {
    mexErrMsgIdAndTxt("cocoProblemGetName:nrhs", "One input required.");
  }
  /* get the problem */
  ref = (size_t *)mxGetData(prhs[0]);
  pb = (coco_problem_t *)(*ref);
  /* call coco_problem_get_name(...) */
  res = coco_problem_get_name(pb);
  /* prepare the return value */
  plhs[0] = mxCreateString(res);
}

void cocoProblemGetNumberOfObjectives(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  size_t *ref;
  coco_problem_t *problem = nullptr;
  const mwSize dims[2] = {1, 1};
  size_t *res;

  /* check for proper number of arguments */
  if (nrhs != 1) {
    mexErrMsgIdAndTxt("cocoProblemGetNumberOfObjectives:nrhs", "One input required.");
  }
  /* get the problem */
  ref = (size_t *)mxGetData(prhs[0]);
  problem = (coco_problem_t *)(*ref);
  /* prepare the return value */
  plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
  res = (size_t *)mxGetData(plhs[0]);
  res[0] = coco_problem_get_number_of_objectives(problem);
}

void cocoProblemGetNumberOfConstraints(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  const mwSize dims[2] = {1, 1};
  size_t *res;

  UNPACK_PROBLEM(0);

  /* prepare the return value */
  plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
  res = (size_t *)mxGetData(plhs[0]);
  res[0] = coco_problem_get_number_of_constraints(problem);
}

void cocoProblemGetNumberOfIntegerVariables(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  const mwSize dims[2] = {1, 1};
  size_t *res;

  UNPACK_PROBLEM(0);

  /* prepare the return value */
  plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
  res = (size_t *)mxGetData(plhs[0]);
  res[0] = coco_problem_get_number_of_integer_variables(problem);
}

void cocoProblemGetSmallestValuesOfInterest(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  size_t nb_dim;
  size_t i;
  const double *res;
  double *v; /* intermediate variable that aloows to set plhs[0] */

  UNPACK_PROBLEM(0);

  nb_dim = coco_problem_get_dimension(problem);
  plhs[0] = mxCreateDoubleMatrix(1, nb_dim, mxREAL);
  v = mxGetPr(plhs[0]);
  /* call coco_problem_get_smallest_values_of_interest(...) */
  res = coco_problem_get_smallest_values_of_interest(problem);
  for (i = 0; i < nb_dim; i++) {
    v[i] = res[i];
  }
}

void cocoProblemGetBestValue(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  size_t *ref;
  coco_problem_t *problem = nullptr;
  double res;

  /* Check for proper number of input and arguments */
  if (nrhs != 1) {
    mexErrMsgIdAndTxt("cocoProblemGetBestValue:nrhs", "One input required.");
  }
  if (nlhs < 1) {
    mexErrMsgIdAndTxt("cocoProblemGetBestValue:nlhs", "One output argument required.");
  }

  /* Get the problem */
  ref = (size_t *)mxGetData(prhs[0]);
  problem = (coco_problem_t *)(*ref);

  /* Return best value objective value */
  res = coco_problem_get_best_value(problem);
  plhs[0] = mxCreateDoubleScalar(res);
}

void cocoProblemIsValid(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* Cannot use UNPACK_PROBLEM because it checks for NULL values. */
  UNPACK_THING(coco_problem_t, problem, 0);

  plhs[0] = mxCreateLogicalScalar(problem != nullptr);
}

void cocoProblemRemoveObserver(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  coco_problem_t *unobservedproblem;

  /* check for proper number of arguments */
  if (nrhs != 2) {
    mexErrMsgIdAndTxt("cocoProblemREmoveObserver:nrhs", "Two inputs required.");
  }

  UNPACK_PROBLEM(0);
  UNPACK_OBSERVER(1);

  /* call coco_problem_remove_observer() */
  unobservedproblem = coco_problem_remove_observer(problem, observer);

  RETURN_POINTER(unobservedproblem, 0);
}

void cocoSetLogLevel(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  char *level;
  const char *res;

  /* check for proper number of arguments */
  if (nrhs != 1) {
    mexErrMsgIdAndTxt("cocoSetLogLevel:nrhs", "One input required.");
  }
  /* get the problem */
  level = mxArrayToString(prhs[0]);
  /* call coco_set_log_level(...) */
  res = coco_set_log_level(level);
  /* prepare the return value */
  plhs[0] = mxCreateString(res);
}

void cocoSuite(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  char *suite_name;
  char *suite_instance;
  char *suite_options;
  coco_suite_t *suite = nullptr;

  /* check for proper number of arguments */
  if (nrhs != 3) {
    mexErrMsgIdAndTxt("cocoSuite:nrhs", "Three inputs required.");
  }
  /* get suite_name */
  suite_name = mxArrayToString(prhs[0]);
  /* get suite_instance */
  suite_instance = mxArrayToString(prhs[1]);
  /* get suite_options */
  suite_options = mxArrayToString(prhs[2]);
  /* call coco_suite() */
  suite = coco_suite(suite_name, suite_instance, suite_options);

  RETURN_POINTER(suite, 0);
}

void cocoSuiteFree(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  UNPACK_SUITE(0);
  coco_suite_free(suite);
}

void cocoSuiteGetNextProblem(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  UNPACK_SUITE(0);
  UNPACK_OBSERVER(1);

  coco_problem_t *problem = coco_suite_get_next_problem(suite, observer);

  RETURN_POINTER(problem, 0);
}

void cocoSuiteGetProblem(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs != 2) {
    mexErrMsgIdAndTxt("cocoSuiteGetProblem:nrhs", "Two inputs required.\n Try \'help cocoSuiteGetProblem.m\'.");
  }
  UNPACK_SUITE(0);
  size_t findex = (size_t)mxGetScalar(prhs[1]);

  coco_problem_t *problem = coco_suite_get_problem(suite, findex);

  RETURN_POINTER(problem, 0);
}

/** @brief The gateway function, calling all Coco functionality from Matlab
 *
 * Called as
 *
 *   out = cocoCall('functionName', arguments).
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  char *cocofunction;
  int i;

  /* check for proper number of arguments */
  if (nrhs < 2) {
    mexErrMsgIdAndTxt("cocoCall:nrhs", "At least two inputs required.");
  }
  /* Get the function string, indicating which function to call. */
  cocofunction = mxArrayToString(prhs[0]);
  /* Convert the given string to lower case for easier case handling
       and more flexibility for the user (code copied from stackoverflow).
  */
  for (i = 0; cocofunction[i]; i++) {
    cocofunction[i] = tolower(cocofunction[i]);
  }

  /* Now the big 'switch' accessing all supported Coco functions. */
  if (strcmp(cocofunction, "cocoevaluatefunction") == 0) {
    cocoEvaluateFunction(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoevaluateconstraint") == 0) {
    cocoEvaluateConstraint(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocorecommendsolution") == 0) {
    cocoRecommendSolution(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoobserver") == 0) {
    cocoObserver(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoobserverfree") == 0) {
    cocoObserverFree(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoobserversignalrestart") == 0) {
    cocoObserverSignalRestart(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemaddobserver") == 0) {
    cocoProblemAddObserver(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemfinaltargethit") == 0) {
    cocoProblemFinalTargetHit(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemfree") == 0) {
    cocoProblemFree(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemgetdimension") == 0) {
    cocoProblemGetDimension(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemgetevaluations") == 0) {
    cocoProblemGetEvaluations(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemgetevaluationsconstraints") == 0) {
    cocoProblemGetEvaluationsConstraints(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemgetid") == 0) {
    cocoProblemGetId(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemgetinitialsolution") == 0) {
    cocoProblemGetInitialSolution(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemgetlargestfvaluesofinterest") == 0) {
    cocoProblemGetLargestFValuesOfInterest(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemgetlargestvaluesofinterest") == 0) {
    cocoProblemGetLargestValuesOfInterest(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemgetname") == 0) {
    cocoProblemGetName(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemgetnumberofobjectives") == 0) {
    cocoProblemGetNumberOfObjectives(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemgetnumberofconstraints") == 0) {
    cocoProblemGetNumberOfConstraints(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemgetnumberofintegervariables") == 0) {
    cocoProblemGetNumberOfIntegerVariables(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemgetsmallestvaluesofinterest") == 0) {
    cocoProblemGetSmallestValuesOfInterest(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemgetbestvalue") == 0) {
    cocoProblemGetBestValue(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemisvalid") == 0) {
    cocoProblemIsValid(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocoproblemremoveobserver") == 0) {
    cocoProblemRemoveObserver(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocosetloglevel") == 0) {
    cocoSetLogLevel(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocosuite") == 0) {
    cocoSuite(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocosuitefree") == 0) {
    cocoSuiteFree(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocosuitegetnextproblem") == 0) {
    cocoSuiteGetNextProblem(nlhs, plhs, nrhs - 1, prhs + 1);
  } else if (strcmp(cocofunction, "cocosuitegetproblem") == 0) {
    cocoSuiteGetProblem(nlhs, plhs, nrhs - 1, prhs + 1);
  } else {
    coco_warning("Function string '%s' not supported", cocofunction);
  }

  mxFree(cocofunction);
}
