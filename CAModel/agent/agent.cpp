#include <cstring>

#include "cam_agent.h"


/*
 * =============================================================================
 * Agent constructor
 * =============================================================================
 */
Agent::Agent (const int x, const int y, const AgentType tp, struct cam_rd_sys* rds)
{
  /* coordinates */
  set_coordinates(x, y);

  /* associate rd_sys */
  rd_sys = rds;

  /* age and lifespan */
  age = 0.0;
  lifespan = 0.0;
  max_lifespan = 0.0;

  /* set agent type */
  set_agent_type(tp);

  /* conversion */
  flag_conversion = false;

  /* bursting */
  flag_bursting = false;

  /* apoptosis */
  is_being_apoptosed = false;
  apoptosis_start_time = 0.0;
  apoptosis_timescale = 0.0;
  num_cells_apoptosing = 0;
  is_apoptosing_epithelial = false;
  apoptosing_agent = nullptr;

  /* phagocytosis */
  is_being_phagocytosed = false;
  phagocytosis_start_time = 0.0;
  phagocytosis_timescale = 0.0;
  is_phagocytosing_epithelial = false;
  phagocytosing_agent = nullptr;
  phagocytosed_by = nullptr;
}


/*
 * =============================================================================
 * Agent destructor
 * =============================================================================
 */
Agent::~Agent ()
{
  if (rd_sys != nullptr) {
    rd_sys = nullptr;
  }
}
