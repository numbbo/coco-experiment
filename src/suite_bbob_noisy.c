/**
 * @file suite_bbob_noisy.c
 * @brief Implementation of the bbob suite containing 24 noiseless single-objective functions in 6
 * dimensions.
 */

#include "suite_bbob_noisy_utilities.c"

#include "f_different_powers.c"
#include "f_ellipsoid.c"
#include "f_gallagher.c"
#include "f_griewank_rosenbrock.c"
#include "f_griewank_rosenbrock.c"
#include "f_rosenbrock.c"
#include "f_schaffers.c"
#include "f_sphere.c"
#include "f_step_ellipsoid.c"

static coco_suite_t *coco_suite_allocate(const char *suite_name,
                                         const size_t number_of_functions,
                                         const size_t number_of_dimensions,
                                         const size_t *dimensions,
                                         const char *default_instances);

/**
 * @brief Sets the dimensions and default instances for the bbob suite.
 */
static coco_suite_t *suite_bbob_noisy_initialize(void) {

  coco_suite_t *suite;
  const size_t dimensions[] = { 2, 3, 5, 10, 20, 40 };
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  /* IMPORTANT: Make sure to change the default instance for every new workshop! */
  suite = coco_suite_allocate("bbob-noisy", 30, num_dimensions, dimensions, "year:2009");

  return suite;
}



/**
 * @brief Sets the instances associated with years for the bbob suite.
 */
static const char *suite_bbob_noisy_get_instances_by_year(const int year) {

  if (year <= 2009) {
    return "1-15";
  }
  else {
    coco_error("suite_bbob_noisy_get_instances_by_year(): year %d not defined for suite bbob-noisy", year);
    return NULL;
  }
}


/**
 * @brief Creates and returns a BBOB problem without needing the actual bbob suite.
 *
 * Useful for other suites as well (see for example suite_biobj.c).
 */
static coco_problem_t *coco_get_bbob_noisy_problem(const size_t function,
                                             const size_t dimension,
                                             const size_t instance) {
  coco_problem_t *problem = NULL;

  const int verbose = 0;
  const char *problem_id_template = "bbob_noisy_f%lu_i%02lu_d%02lu";
  const char *problem_name_template = "BBOB-NOISY suite problem f%lu instance %lu in %luD";

  long rseed;
  const long rseed_1  = (long) (1 + 10000 * instance);
  const long rseed_8  = (long) (8 + 10000 * instance);
  const long rseed_7  = (long) (7 + 10000 * instance);
  const long rseed_10 = (long) (10 + 10000 * instance);
  const long rseed_14 = (long) (14 + 10000 * instance);
  const long rseed_17 = (long) (17 + 10000 * instance);
  const long rseed_19 = (long) (19 + 10000 * instance);
  const long rseed_21 = (long) (21 + 10000 * instance);

  const uint32_t coco_seed = (uint32_t) (31 + 1000 + instance);

  const size_t function_idx = function + 100;

  if (function_idx == 101){
    double *distribution_theta = coco_allocate_vector((const size_t) 1);
    distribution_theta[0] = 0.01;
    rseed = rseed_1;
    problem = coco_problem_allocate_bbob_wrap_noisy_gaussian(
        f_sphere_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  102){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    double alpha, beta;
    alpha = 0.01 * (0.49 + 1 / (double) dimension);
    beta = 0.01;
    distribution_theta[0] = alpha;
    distribution_theta[1] = beta;
    rseed = rseed_1;
    problem = coco_problem_allocate_bbob_wrap_noisy_uniform(
        f_sphere_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  103){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    distribution_theta[0] = 0.01;
    distribution_theta[1] = 0.05;
    rseed = rseed_1;
    problem = coco_problem_allocate_bbob_wrap_noisy_cauchy(
        f_sphere_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  104){
    double *distribution_theta = coco_allocate_vector((const size_t) 1);
    distribution_theta[0] = 0.01; 
    rseed = rseed_8;
    problem = coco_problem_allocate_bbob_wrap_noisy_gaussian(
        f_rosenbrock_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  105){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    double alpha, beta;
    alpha = 0.01 * (0.49 + 1 / (double) dimension);
    beta = 0.01;
    distribution_theta[0] = alpha;
    distribution_theta[1] = beta;
    rseed = rseed_8;
    problem = coco_problem_allocate_bbob_wrap_noisy_uniform(
        f_rosenbrock_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  106){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    distribution_theta[0] = 0.01; 
    distribution_theta[1] = 0.05; 
    rseed = rseed_8;
    problem = coco_problem_allocate_bbob_wrap_noisy_cauchy(
        f_rosenbrock_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  107){
    double *distribution_theta = coco_allocate_vector((const size_t) 1);
    distribution_theta[0] = 1.0; 
    rseed = rseed_1;
    problem = coco_problem_allocate_bbob_wrap_noisy_gaussian(
        f_sphere_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  108){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    double alpha, beta;
    alpha = 0.49 + 1.0 / (double) dimension;
    beta = 1.0;
    distribution_theta[0] = alpha;
    distribution_theta[1] = beta;
    rseed = rseed_1;
    problem = coco_problem_allocate_bbob_wrap_noisy_uniform(
        f_sphere_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  109){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    distribution_theta[0] = 1.0;
    distribution_theta[1] = 0.2;
    rseed = rseed_1;
    problem = coco_problem_allocate_bbob_wrap_noisy_cauchy(
        f_sphere_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  110){
    double *distribution_theta = coco_allocate_vector((const size_t) 1);
    distribution_theta[0] = 1.0; 
    rseed = rseed_8;
    problem = coco_problem_allocate_bbob_wrap_noisy_gaussian(
        f_rosenbrock_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  111){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    double alpha, beta;
    alpha = 0.49 + 1.0 / (double) dimension;
    beta = 1.0;
    distribution_theta[0] = alpha;
    distribution_theta[1] = beta;
    rseed = rseed_8;
    problem = coco_problem_allocate_bbob_wrap_noisy_uniform(
        f_rosenbrock_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  112){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    distribution_theta[0]= 1.0;
    distribution_theta[1]= 0.2;
    rseed = rseed_8;
    problem = coco_problem_allocate_bbob_wrap_noisy_cauchy(
        f_rosenbrock_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  113){
    double *distribution_theta = coco_allocate_vector((const size_t) 1);
    distribution_theta[0] = 1.0; 
    rseed = rseed_7;
    problem = coco_problem_allocate_bbob_wrap_noisy_gaussian(
        f_step_ellipsoid_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  114){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    double alpha, beta;
    alpha = 0.49 + 1.0 / (double) dimension;
    beta = 1.0;
    distribution_theta[0] = alpha;
    distribution_theta[1] = beta;
    rseed = rseed_7;
    problem = coco_problem_allocate_bbob_wrap_noisy_uniform(
        f_step_ellipsoid_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  115){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    distribution_theta[0] = 1.0; 
    distribution_theta[1]= 0.2; 
    rseed = rseed_7;
    problem = coco_problem_allocate_bbob_wrap_noisy_cauchy(
        f_step_ellipsoid_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  116){
    double *distribution_theta = coco_allocate_vector((const size_t) 1);
    distribution_theta[0] = 1.0; 
    rseed = rseed_10;
    double condition = 1.0e4;
    problem = coco_problem_allocate_bbob_wrap_noisy_gaussian_conditioned(
        f_ellipsoid_rotated_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        condition,
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  117){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    double alpha, beta;
    alpha = 0.49 + 1.0 / (double) dimension;
    beta = 1.0;
    distribution_theta[0] = alpha; 
    distribution_theta[1] = beta; 
    rseed = rseed_10;
    double condition = 1.0e4;
    problem = coco_problem_allocate_bbob_wrap_noisy_uniform_conditioned(
        f_ellipsoid_rotated_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        condition,
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  118){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    distribution_theta[0]= 1.0; 
    distribution_theta[1]= 0.2; 
    rseed = rseed_10;
    double condition = 1.0e4;
    problem = coco_problem_allocate_bbob_wrap_noisy_cauchy_conditioned(
        f_ellipsoid_rotated_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed,
        condition, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  119){
    double *distribution_theta = coco_allocate_vector((const size_t) 1);
    distribution_theta[0] = 1.0; 
    rseed = rseed_14;
    problem = coco_problem_allocate_bbob_wrap_noisy_gaussian(
        f_different_powers_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  120){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    double alpha, beta;
    alpha = 0.49 + 1.0 / (double) dimension;
    beta = 1.0;
    distribution_theta[0] = alpha; 
    distribution_theta[1] = beta; 
    rseed = rseed_14;
    problem = coco_problem_allocate_bbob_wrap_noisy_uniform(
        f_different_powers_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  121){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    distribution_theta[0]= 1.0; 
    distribution_theta[1]= 0.2; 
    rseed = rseed_14;
    problem = coco_problem_allocate_bbob_wrap_noisy_cauchy(
        f_different_powers_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  122){
    double *distribution_theta = coco_allocate_vector((const size_t) 1);
    distribution_theta[0]= 1.0; 
    double conditioning = 10; 
    rseed = rseed_17;
    problem = coco_problem_allocate_bbob_wrap_noisy_gaussian_conditioned(
        f_schaffers_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed,
        conditioning, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  123){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    double alpha, beta;
    alpha = 0.49 + 1.0 / (double) dimension;
    beta = 1.0;
    distribution_theta[0] = alpha; 
    distribution_theta[1] = beta; 
    double conditioning = 10;
    rseed = rseed_17;
    problem = coco_problem_allocate_bbob_wrap_noisy_uniform_conditioned(
        f_schaffers_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed,
        conditioning, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  124){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    distribution_theta[0] = 1.0; 
    distribution_theta[1] = 0.2; 
    double conditioning = 10; 
    rseed = rseed_17;
    problem = coco_problem_allocate_bbob_wrap_noisy_cauchy_conditioned(
        f_schaffers_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed_17,
        conditioning, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  125){
    double *distribution_theta = coco_allocate_vector((const size_t) 1);
    distribution_theta[0]= 1.0; 
    rseed = rseed_19;
    double facftrue = 1.0;
    problem = coco_problem_allocate_bbob_wrap_noisy_gaussian_conditioned(
        f_griewank_rosenbrock_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed,
        facftrue, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  126){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    double alpha, beta;
    alpha = 0.49 + 1.0 / (double) dimension;
    beta = 1.0;
    distribution_theta[0] = alpha; 
    distribution_theta[1] = beta; 
    rseed = rseed_19;
    double facftrue = 1.0;
    problem = coco_problem_allocate_bbob_wrap_noisy_uniform_conditioned(
        f_griewank_rosenbrock_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed,
        facftrue, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  127){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    distribution_theta[0] = 1.0; 
    distribution_theta[1] = 0.2; 
    rseed = rseed_19;
    double facftrue = 1.0;
    problem = coco_problem_allocate_bbob_wrap_noisy_cauchy_conditioned(
        f_griewank_rosenbrock_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed, 
        facftrue, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  128){
    double *distribution_theta = coco_allocate_vector((const size_t) 1);
    distribution_theta[0]= 1.0;
    size_t n_peaks = 101; 
    rseed = rseed_21;
    problem = coco_problem_allocate_bbob_wrap_noisy_gaussian_gallagher(
        f_gallagher_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed,
        n_peaks,
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  129){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    double alpha, beta;
    alpha = 0.49 + 1.0 / (double) dimension;
    beta = 1.0;
    distribution_theta[0] = alpha; 
    distribution_theta[1] = beta; 
    size_t n_peaks = 101;
    rseed = rseed_21;
    problem = coco_problem_allocate_bbob_wrap_noisy_uniform_gallagher(
        f_gallagher_bbob_problem_allocate, 
        function_idx, 
        dimension, 
        instance, 
        rseed,
        n_peaks, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else if (function_idx ==  130){
    double *distribution_theta = coco_allocate_vector((const size_t) 2);
    distribution_theta[0] = 1.0; 
    distribution_theta[1] = 0.2; 
    size_t n_peaks = 101; 
    rseed = rseed_21;
    problem = coco_problem_allocate_bbob_wrap_noisy_cauchy_gallagher(
        f_gallagher_bbob_problem_allocate, 
        function_idx,
        dimension, 
        instance, 
        rseed,
        n_peaks, 
        problem_id_template, 
        problem_name_template,
        coco_seed, 
        distribution_theta
    );
  } else {
    coco_error("coco_get_bbob_noisy_problem(): cannot retrieve problem f%lu instance %lu in %luD",
    		(unsigned long) function_idx, (unsigned long) instance, (unsigned long) dimension);
    return NULL; /* Never reached */
  }
  if (verbose == 1){    
    printf("\nf%lu_i%02lu_d%02lu, optimal_f_value: %f \n", function_idx, instance, dimension , problem -> best_value[0]);
    printf("x: {\n");
    for(size_t dim = 0; dim < problem -> number_of_variables - 1; dim ++){
        printf("\t%f, \n", *(problem -> best_parameter + dim));
        
    }
    printf("\t%f\n", *(problem -> best_parameter + problem -> number_of_variables - 1));
    printf("}");
  }
  if (verbose == 2) {
    printf("f%lu_i%02lu_d%02lu, rseed: %lu \n", function_idx, instance, dimension , rseed);

  }
  return problem;
}

/**
 * @brief Returns the problem from the bbob suite that corresponds to the given parameters.
 *
 * @param suite The COCO suite.
 * @param function_idx Index of the function (starting from 0).
 * @param dimension_idx Index of the dimension (starting from 0).
 * @param instance_idx Index of the instance (starting from 0).
 * @return The problem that corresponds to the given parameters.
 */
static coco_problem_t *suite_bbob_noisy_get_problem( coco_suite_t *suite,
                                              const size_t function_idx,
                                              const size_t dimension_idx,
                                              const size_t instance_idx) {

  coco_problem_t *problem = NULL;

  const size_t function = suite->functions[function_idx];
  const size_t dimension = suite->dimensions[dimension_idx];
  const size_t instance = suite->instances[instance_idx];

  problem = coco_get_bbob_noisy_problem(function, dimension, instance);

  problem -> suite_dep_function = function;
  problem -> suite_dep_instance = instance;
  problem -> suite_dep_index = coco_suite_encode_problem_index(suite, function_idx, dimension_idx, instance_idx);
  return problem;
}
