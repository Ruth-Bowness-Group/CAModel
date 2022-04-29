#include <algorithm>

#include "cam_util.h"
#include "cam_parameters.h"


/*
 * =============================================================================
 * Get the moore neighbourhood around a given coordinate
 * =============================================================================
 */

/* moore neighbourhood around coords xy */
std::vector<std::pair<int, int>>
cam_util_neighbourhood_moore (const std::pair<int,int> xy, const int level)
{
  int i = 0, j = 0;
  std::vector<std::pair<int, int>> neighbours;

  int x = xy.first;
  int y = xy.second;

  for (i = -level; i <= level; i++) {
    for (j = -level; j <= level; j++) {
      /* Check that neighbour is not the original loc and is on the grid */
      if ( !(i == 0 && j == 0) && ((x + i) >= 0) && ((x + i) < grid_size) &&
           ((y + j) >= 0) && ((y + j) < grid_size) ) {
        /* Create a pair of coordinates and add it to the vector */
        std::pair<int, int> neighbour = std::make_pair(x + i, y + j);
        neighbours.push_back(neighbour);
      }
    }
  }

  /* Shuffle the neighbours */
  std::shuffle(neighbours.begin(), neighbours.end(), rng);

  return neighbours;
}


/* moore neighbourhood around coords (x,y) */
std::vector<std::pair<int, int>>
cam_util_neighbourhood_moore (const int x, const int y, const int level)
{
  int i = 0, j = 0;
  std::vector<std::pair<int, int>> neighbours;

  for (i = -level; i <= level; i++) {
    for (j = -level; j <= level; j++) {
      /* Check that neighbour is not the original loc and is on the grid */
      if ( !(i == 0 && j == 0) && ((x + i) >= 0) && ((x + i) < grid_size) &&
           ((y + j) >= 0) && ((y + j) < grid_size) ) {
        /* Create a pair of coordinates and add it to the vector */
        std::pair<int, int> neighbour = std::make_pair(x + i, y + j);
        neighbours.push_back(neighbour);
      }
    }
  }

  /* Shuffle the neighbours */
  std::shuffle(neighbours.begin(), neighbours.end(), rng);

  return neighbours;
}


double
cam_util_tanh_function (const double val, const double shift, const double scale, const bool too_large)
{
  if (too_large)
    return 1.0;

  double arg = scale * (val - shift);
  return (0.5 * (1.0 + tanh(arg)));
}


double
cam_util_tanh_decay_function (const double val, const double shift, const double scale)
{
  double arg = scale * (val - shift);
  return (0.5 * (1.0 - tanh(arg)));
}
