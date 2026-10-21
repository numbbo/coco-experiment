/**
 * @file suite_bbob_problems.cpp
 * @brief The problems of the bbob, bbob-boxed and bbob-largescale suites.
 */

#include "coco.h"

#include "f_attractive_sector.cpp"
#include "f_bent_cigar.cpp"
#include "f_bent_cigar_generalized.cpp"
#include "f_bueche_rastrigin.cpp"
#include "f_different_powers.cpp"
#include "f_discus.cpp"
#include "f_discus_generalized.cpp"
#include "f_ellipsoid.cpp"
#include "f_gallagher.cpp"
#include "f_griewank_rosenbrock.cpp"
#include "f_katsuura.cpp"
#include "f_linear_slope.cpp"
#include "f_lunacek_bi_rastrigin.cpp"
#include "f_rastrigin.cpp"
#include "f_rosenbrock.cpp"
#include "f_schaffers.cpp"
#include "f_schwefel.cpp"
#include "f_sharp_ridge.cpp"
#include "f_sharp_ridge_generalized.cpp"
#include "f_sphere.cpp"
#include "f_step_ellipsoid.cpp"
#include "f_weierstrass.cpp"

namespace coco {

/**
 * @brief Creates and returns the problem with the given function, dimension and instance of a bbob suite variant.
 *
 * The large-scale variant uses permuted block rotations instead of full rotations for the non-separable functions.
 * The id and name templates are formatted with the function, instance and dimension.
 */
static ProblemPtr MakeBbobProblem(const size_t function, const size_t dimension, const size_t instance,
                                  const BbobVariant variant, const char *problem_id_template,
                                  const char *problem_name_template) {
  const long rseed = static_cast<long>(function + 10000 * instance);
  const long rseed_3 = static_cast<long>(3 + 10000 * instance);
  const long rseed_17 = static_cast<long>(17 + 10000 * instance);
  const bool largescale = variant == BbobVariant::LARGESCALE;

  switch (function) {
  case 1:
    return std::make_unique<BbobSphere>(function, dimension, instance, rseed, variant, problem_id_template,
                                        problem_name_template);
  case 2:
    return std::make_unique<BbobEllipsoid>(function, dimension, instance, rseed, variant, problem_id_template,
                                           problem_name_template);
  case 3:
    return std::make_unique<BbobRastrigin>(function, dimension, instance, rseed, variant, problem_id_template,
                                           problem_name_template);
  case 4:
    return std::make_unique<BbobBuecheRastrigin>(function, dimension, instance, rseed_3, variant, problem_id_template,
                                                 problem_name_template);
  case 5:
    return std::make_unique<BbobLinearSlope>(function, dimension, instance, rseed, variant, problem_id_template,
                                             problem_name_template);
  case 6:
    if (largescale)
      return std::make_unique<BbobAttractiveSectorPermblockdiag>(function, dimension, instance, rseed, variant,
                                                                 problem_id_template, problem_name_template);
    return std::make_unique<BbobAttractiveSector>(function, dimension, instance, rseed, variant, problem_id_template,
                                                  problem_name_template);
  case 7:
    if (largescale)
      return std::make_unique<BbobStepEllipsoidPermblockdiag>(function, dimension, instance, rseed, variant,
                                                              problem_id_template, problem_name_template);
    return std::make_unique<BbobStepEllipsoid>(function, dimension, instance, rseed, variant, 1.0, problem_id_template,
                                               problem_name_template);
  case 8:
    return std::make_unique<BbobRosenbrock>(function, dimension, instance, rseed, variant, problem_id_template,
                                            problem_name_template);
  case 9:
    if (largescale)
      return std::make_unique<BbobRosenbrockPermblockdiag>(function, dimension, instance, rseed, variant,
                                                           problem_id_template, problem_name_template);
    return std::make_unique<BbobRosenbrockRotated>(function, dimension, instance, rseed, problem_id_template,
                                                   problem_name_template);
  case 10:
    if (largescale)
      return std::make_unique<BbobEllipsoidPermblockdiag>(function, dimension, instance, rseed, variant,
                                                          problem_id_template, problem_name_template);
    return std::make_unique<BbobEllipsoidRotated>(function, dimension, instance, rseed, variant, 1.0e6,
                                                  problem_id_template, problem_name_template);
  case 11:
    if (largescale)
      return std::make_unique<BbobDiscusGeneralizedPermblockdiag>(function, dimension, instance, rseed, variant,
                                                                  problem_id_template, problem_name_template);
    return std::make_unique<BbobDiscus>(function, dimension, instance, rseed, variant, problem_id_template,
                                        problem_name_template);
  case 12:
    if (largescale)
      return std::make_unique<BbobBentCigarGeneralizedPermblockdiag>(function, dimension, instance, rseed, variant,
                                                                     problem_id_template, problem_name_template);
    return std::make_unique<BbobBentCigar>(function, dimension, instance, rseed, variant, problem_id_template,
                                           problem_name_template);
  case 13:
    if (largescale)
      return std::make_unique<BbobSharpRidgeGeneralizedPermblockdiag>(function, dimension, instance, rseed, variant,
                                                                      problem_id_template, problem_name_template);
    return std::make_unique<BbobSharpRidge>(function, dimension, instance, rseed, variant, problem_id_template,
                                            problem_name_template);
  case 14:
    if (largescale)
      return std::make_unique<BbobDifferentPowersPermblockdiag>(function, dimension, instance, rseed, variant,
                                                                problem_id_template, problem_name_template);
    return std::make_unique<BbobDifferentPowers>(function, dimension, instance, rseed, variant, problem_id_template,
                                                 problem_name_template);
  case 15:
    if (largescale)
      return std::make_unique<BbobRastriginRotatedPermblockdiag>(function, dimension, instance, rseed, variant,
                                                                 problem_id_template, problem_name_template);
    return std::make_unique<BbobRastriginRotated>(function, dimension, instance, rseed, variant, problem_id_template,
                                                  problem_name_template);
  case 16:
    if (largescale)
      return std::make_unique<BbobWeierstrassPermblockdiag>(function, dimension, instance, rseed, variant,
                                                            problem_id_template, problem_name_template);
    return std::make_unique<BbobWeierstrass>(function, dimension, instance, rseed, variant, problem_id_template,
                                             problem_name_template);
  case 17:
    if (largescale)
      return std::make_unique<BbobSchaffersPermblockdiag>(function, dimension, instance, rseed, variant, 10.0,
                                                          problem_id_template, problem_name_template);
    return std::make_unique<BbobSchaffers>(function, dimension, instance, rseed, variant, 10.0, 10.0,
                                           problem_id_template, problem_name_template);
  case 18:
    if (largescale)
      return std::make_unique<BbobSchaffersPermblockdiag>(function, dimension, instance, rseed_17, variant, 1000.0,
                                                          problem_id_template, problem_name_template);
    return std::make_unique<BbobSchaffers>(function, dimension, instance, rseed_17, variant, 1000.0, 10.0,
                                           problem_id_template, problem_name_template);
  case 19:
    if (largescale)
      return std::make_unique<BbobGriewankRosenbrockPermblockdiag>(function, dimension, instance, rseed, 10.0,
                                                                   problem_id_template, problem_name_template);
    return std::make_unique<BbobGriewankRosenbrock>(function, dimension, instance, rseed, 10.0, problem_id_template,
                                                    problem_name_template);
  case 20:
    return std::make_unique<BbobSchwefel>(function, dimension, instance, rseed, variant, problem_id_template,
                                          problem_name_template);
  case 21:
    if (largescale)
      return std::make_unique<BbobGallagherPermblockdiag>(function, dimension, instance, rseed, 101,
                                                          problem_id_template, problem_name_template);
    return std::make_unique<BbobGallagher>(function, dimension, instance, rseed, variant, static_cast<size_t>(101), 1.0,
                                           problem_id_template, problem_name_template);
  case 22:
    if (largescale)
      return std::make_unique<BbobGallagherPermblockdiag>(function, dimension, instance, rseed, 21, problem_id_template,
                                                          problem_name_template);
    return std::make_unique<BbobGallagher>(function, dimension, instance, rseed, variant, static_cast<size_t>(21), 1.0,
                                           problem_id_template, problem_name_template);
  case 23:
    if (largescale)
      return std::make_unique<BbobKatsuuraPermblockdiag>(function, dimension, instance, rseed, variant,
                                                         problem_id_template, problem_name_template);
    return std::make_unique<BbobKatsuura>(function, dimension, instance, rseed, variant, problem_id_template,
                                          problem_name_template);
  case 24:
    if (largescale)
      return std::make_unique<BbobLunacekBiRastriginPermblockdiag>(function, dimension, instance, rseed,
                                                                   problem_id_template, problem_name_template);
    return std::make_unique<BbobLunacekBiRastrigin>(function, dimension, instance, rseed, variant, problem_id_template,
                                                    problem_name_template);
  default:
    coco_error("MakeBbobProblem(): cannot retrieve problem f%zu instance %zu in %zuD", function, instance, dimension);
    return nullptr; /* Never reached */
  }
}

} // namespace coco
