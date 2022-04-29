#ifndef __CAM_VESSEL_H__
#define __CAM_VESSEL_H__


#include <vector>

#include "../reactdiff/cam_reactdiff.h"


/* Set the Vessel specific parameters */
int cam_vessel_set_params (YAML::Node* cfg);


/* enum class of Vessel Type */
enum class VesselType
{
  null = -1, // no type
  blood = 0, // blood vessel
  lymph = 1 // lymph vessel
};


class Vessel
{
private:
  /* coordinates */
  int x;
  int y;

  /* area */
  double area;
  double max_area;

  /* vessel type */
  VesselType type;

  /* reaction diffusion system */
  struct cam_rd_sys* rd_sys;

  /* Secrete */

  /* Uptake */
  double uptake_ifn_1;
  double uptake_vir_inf;

  /* set the uptake */
  void set_uptake (const VesselType tp);

public:
  Vessel (const int xc, const int yc, const VesselType tp, struct cam_rd_sys* rds, int* status);
  ~Vessel ();

  /* neighbours */
  std::vector<std::pair<int,int>> immediate_moore;
  std::vector<std::pair<int,int>> neighbours;

  /* coordinates */
  std::pair<int,int> get_coordinates ();
  void set_coordinates (const int cx, const int cy);

  /* get area */
  double get_max_area ();
  bool check_area (const double area_to_add);

  /* vessel type */
  VesselType get_type ();
  void set_type (const VesselType tp);

  /* output */
  int get_code ();
};


/*
 * =============================================================================
 * Inline functions
 * =============================================================================
 */

/* get the coordinates */
inline std::pair<int,int>
Vessel::get_coordinates ()
{
  return std::make_pair(x, y);
}

/* set the coordinates */
inline void
Vessel::set_coordinates (const int cx, const int cy)
{
  x = cx;
  y = cy;
}

/* check whether agent can be added based on area */
inline bool
Vessel::check_area (const double area_to_add)
{
  double tmp_area = area + area_to_add;
  return ( (CAM_CMP(tmp_area, max_area, CAM_TOL)) ? true : false ); /* area + area_to_add < max_area */
}

/* get the area */
inline double
Vessel::get_max_area ()
{
  return max_area;
}

/* get vessel type */
inline VesselType
Vessel::get_type ()
{
  return type;
}

/* set the vessel type */
inline void
Vessel::set_type (const VesselType tp)
{
  type = tp;
}

/* get the output grid code */
inline int
Vessel::get_code ()
{
  return static_cast<int>(type);
}

#endif /* __CAM_VESSEL_H__ */
