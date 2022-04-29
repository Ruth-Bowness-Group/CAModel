#ifndef __CAM_RANDOM_H__
#define __CAM_RANDOM_H__

#include <random>

#include "cam_parameters.h"


/*
 * Random number generator.
 * Uniform distributions. Call with rng to generate random number in range.
 *
 * For example, to assign a random number from [0,1] to 'a' use:
 *    double a = distribution_0_1(rng);
 */
static int seed = 1580155444;
static auto rng = std::default_random_engine(seed);

static std::uniform_real_distribution<double> distribution_0_1(0.0, 1.0);
static std::uniform_real_distribution<double> distribution_neg1_1(-1.0, 1.0);
static std::uniform_int_distribution<int> distribution_0_grid_size(0, grid_size-1);

#endif /* __CAM_RANDOM_H__ */
