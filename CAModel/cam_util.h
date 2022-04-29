#ifndef __CAM_UTIL_H__
#define __CAM_UTIL_H__

#include <vector>
#include <random>

#include "cam_random.h"


/* moore neighbourhood */
std::vector<std::pair<int, int>> cam_util_neighbourhood_moore (const std::pair<int,int> xy, const int level);
std::vector<std::pair<int, int>> cam_util_neighbourhood_moore (const int x, const int y, const int level);

/* hill function */
inline double
cam_util_hill_function (const double val, const double ec50, const int coeff)
{
  if (val < CAM_TOL) {
    return 0.0;
  }

  double numer = pow(val, coeff);
  double denom = pow(ec50, coeff) + numer;

  return (numer / denom);
}

/* tanh functions */
double cam_util_tanh_function (const double val, const double shift, const double scale, const bool too_large);
double cam_util_tanh_decay_function (const double val, const double shift, const double scale);

inline int
cam_util_choose_random_index (const int lo, const int hi)
{
  std::uniform_int_distribution<int> dist_ind(lo, hi);
  int sp_i = dist_ind(rng);

  return sp_i;
}

#endif /* __CAM_UTIL_H__ */
