#ifndef __CAM_LOCATION_H__
#define __CAM_LOCATION_H__

#include "vessel/cam_vessel.h"
#include "epithelial/cam_epithelial.h"

/*
 * =============================================================================
 * Location class
 * =============================================================================
 */
class Location
{
private:
  int x;
  int y;

  bool is_a_vessel;
  Vessel* vessel;

  bool is_a_epithelial;
  Epithelial* epithelial;

public:
  Location (const int xc, const int yc);
  ~Location ();

  /* coordinates */
  std::pair<int,int> get_coordinates ();
  void set_coordinates (const int cx, const int cy);

  /* Vessel */
  void set_vessel (Vessel* ves);
  Vessel* get_vessel ();
  bool is_vessel ();

  /* Epithelial */
  void set_epithelial (Epithelial* epi);
  Epithelial* get_epithelial ();
  bool is_epithelial ();
};


/*
 * =============================================================================
 * Location specific inline functions
 * =============================================================================
 */

/* get the coordinates */
inline std::pair<int,int>
Location::get_coordinates ()
{
  return std::make_pair(x, y);
}

/* set the coordinates */
inline void
Location::set_coordinates (const int cx, const int cy)
{
  x = cx;
  y = cy;
}

/* set the location as a vessel */
inline void
Location::set_vessel (Vessel* ves)
{
  vessel = ves;
  is_a_vessel = ( (ves == nullptr) ? false : true );
  is_a_epithelial = false;
}

/* get the vessel at the location */
inline Vessel*
Location::get_vessel ()
{
  return vessel;
}

/* check whether grid location is a vessel */
inline bool
Location::is_vessel ()
{
  return is_a_vessel;
}

/* set the location as an epithelial cell */
inline void
Location::set_epithelial (Epithelial* epi)
{
  epithelial = epi;
  is_a_epithelial = ( (epi == nullptr) ? false : true );
  is_a_vessel = false;
}

/* get the epithelial at the location */
inline Epithelial*
Location::get_epithelial ()
{
  return epithelial;
}

/* check whether grid location is an epithelial cell */
inline bool
Location::is_epithelial ()
{
  return is_a_epithelial;
}

#endif /* __CAM_LOCATION_H__ */
