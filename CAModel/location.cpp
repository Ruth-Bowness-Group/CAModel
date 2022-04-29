#include "cam_location.h"


/*
 * =============================================================================
 * Location constructor
 * =============================================================================
 */
Location::Location (const int xc, const int yc)
{
  /* set the coordinates */
  set_coordinates(xc, yc);

  /* initialise the vessel and epithelial parts */
  set_vessel(nullptr);
  set_epithelial(nullptr);
}


/*
 * =============================================================================
 * Location destructor
 * =============================================================================
 */
Location::~Location ()
{
  if (vessel != nullptr) {
    vessel = nullptr;
  }

  if (epithelial != nullptr) {
    epithelial = nullptr;
  }
}
