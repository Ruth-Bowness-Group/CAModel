#include "cam_vessel.h"
#include "vessel.h"
#include "../cam_error.h"


/* struct containing vessel specific parameters */
static vessel_params params;


/*
 * =============================================================================
 * Set the vessel specific parameters
 * =============================================================================
 */
int
cam_vessel_set_params (YAML::Node* cfg)
{
  try {
    /* vessel radii */
    params.vessel_radius = (*cfg)["vessel_radius"].as<double>();

    /* cytokine uptake by vessels */
    params.vessel_blood_uptake_ifn_1 = (*cfg)["vessel_blood_uptake_ifn_1"].as<double>();
    params.vessel_lymph_uptake_ifn_1 = (*cfg)["vessel_lymph_uptake_ifn_1"].as<double>();

    /* virion uptake by vessels */
    params.vessel_blood_uptake_vir_inf = (*cfg)["vessel_blood_uptake_vir_inf"].as<double>();
    params.vessel_lymph_uptake_vir_inf = (*cfg)["vessel_lymph_uptake_vir_inf"].as<double>();

  } catch (std::exception &e) {
    CAM_ERROR(e.what(), CAM_ERROR_IO);
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Vessel constructor
 * =============================================================================
 */
Vessel::Vessel (const int xc, const int yc, const VesselType tp, struct cam_rd_sys* rds, int* status)
{
  /* set the coordinates */
  set_coordinates(xc, yc);

  /* associate rd_sys */
  rd_sys = rds;

  /* area */
  area = 0.0;
  max_area = CAM_PI * params.vessel_radius * params.vessel_radius;

  /* set the vessel type */
  set_type(tp);

  /* set the secretion */

  /* set the uptake */
  set_uptake(tp);

  /* update the reaction diffusion uptakes for this epithelial */
  rd_sys->ifn_1->inc_uptake(xc, yc, uptake_ifn_1, true);
  rd_sys->virus_inf->inc_uptake(xc, yc, uptake_vir_inf, true);

  /* successful */
  (*status) = CAM_SUCCESS;
}


/*
 * =============================================================================
 * Vessel destructor
 * =============================================================================
 */
Vessel::~Vessel ()
{
  if (rd_sys != nullptr) {
    rd_sys = nullptr;
  }
}


/*
 * =============================================================================
 * Set the uptake based on the type of vessel
 * =============================================================================
 */
void
Vessel::set_uptake (const VesselType tp)
{
  switch (tp) {
    case VesselType::blood:
      uptake_ifn_1 = params.vessel_blood_uptake_ifn_1;
      uptake_vir_inf = params.vessel_blood_uptake_vir_inf;
      break;
    case VesselType::lymph:
      uptake_ifn_1 = params.vessel_lymph_uptake_ifn_1;
      uptake_vir_inf = params.vessel_lymph_uptake_vir_inf;
      break;
    default:
      uptake_ifn_1 = 0.0;
      uptake_vir_inf = 0.0;
      break;
  }
}
