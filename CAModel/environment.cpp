#include "cam_environment.h"
#include "environment.h"
#include "cam_error.h"
#include "cam_random.h"
#include "cam_util.h"

#include "agent/innate/cam_macrophage.h"
#include "agent/innate/cam_nk.h"
#include "agent/innate/cam_neutrophil.h"
#include "agent/innate/cam_monocyte.h"

/* struct containing environment specific parameters */
static struct env_params params;


/*
 * =============================================================================
 * Set the environment specific parameters
 * =============================================================================
 */
int
cam_env_set_params (YAML::Node* cfg)
{
  int status = CAM_ERROR_MISC;

  try {
    /* Blood vessels */
    params.env_fixed_blood_vessel_placement = (*cfg)["env_fixed_blood_vessel_placement"].as<bool>();
    params.env_initial_number_blood_vessels = (*cfg)["env_initial_number_blood_vessels"].as<int>();
  } catch (std::exception &e) {
    CAM_ERROR("blood vessels", CAM_ERROR_IO);
  }

  try {
    /* Lymphatic vessels */
    params.env_fixed_lymph_vessel_placement = (*cfg)["env_fixed_lymph_vessel_placement"].as<bool>();
    params.env_initial_number_lymph_vessels = (*cfg)["env_initial_number_lymph_vessels"].as<int>();
  } catch (std::exception &e) {
    CAM_ERROR("lymph vessels", CAM_ERROR_IO);
  }

  try {
    /* Virions */
    params.env_fixed_virion_placement = (*cfg)["env_fixed_virion_placement"].as<bool>();
    params.env_initial_number_virions = (*cfg)["env_initial_number_virions"].as<int>();
    params.env_initial_moi = (*cfg)["env_initial_moi"].as<double>();
  } catch (std::exception &e) {
    CAM_ERROR("virions", CAM_ERROR_IO);
  }

  try {
    /* Macrophages */
    params.env_fixed_macrophage_placement = (*cfg)["env_fixed_macrophage_placement"].as<bool>();
    params.env_initial_number_macrophages = (*cfg)["env_initial_number_macrophages"].as<int>();
  } catch (std::exception &e) {
    CAM_ERROR("macs", CAM_ERROR_IO);
  }

  try {
    /* NKs */
    params.env_fixed_nk_placement = (*cfg)["env_fixed_nk_placement"].as<bool>();
    params.env_initial_number_nks = (*cfg)["env_initial_number_nks"].as<int>();
  } catch (std::exception &e) {
    CAM_ERROR("NKs", CAM_ERROR_IO);
  }

  try {
    /* Neutrohils */
    params.env_fixed_neutrophil_placement = (*cfg)["env_fixed_neutrophil_placement"].as<bool>();
    params.env_initial_number_neutrophils = (*cfg)["env_initial_number_neutrophils"].as<int>();
  } catch (std::exception &e) {
    CAM_ERROR("neuts", CAM_ERROR_IO);
  }

  try {
    /* Monocytes */
    params.env_fixed_monocyte_placement = (*cfg)["env_fixed_monocyte_placement"].as<bool>();
    params.env_initial_number_monocytes = (*cfg)["env_initial_number_monocytes"].as<int>();
  } catch (std::exception &e) {
    CAM_ERROR("monos", CAM_ERROR_IO);
  }

  try {
    /* neighbourhood */
    params.env_neighbourhood_max_level = (*cfg)["env_neighbourhood_max_level"].as<int>();
  } catch (std::exception &e) {
    CAM_ERROR("neighbourhoods", CAM_ERROR_IO);
  }

  /* set the vessel specific parameters */
  status = cam_vessel_set_params(cfg);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* set the epithelial specific parameters */
  status = cam_epithelial_set_params(cfg);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* set the reaction diffusion specific parameters */
  status = cam_rd_set_params(cfg);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* set the macrophage specific parameters */
  status = cam_macrophage_set_params(cfg);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* set the NK specific parameters */
  status = cam_nk_set_params(cfg);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* set the neutrophil specific parameters */
  status = cam_neutrophil_set_params(cfg);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* set the monocyte specific parameters */
  status = cam_monocyte_set_params(cfg);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Environment constructor
 * =============================================================================
 */
Environment::Environment (int* status)
{
  int x = 0, y = 0;

  /* counter for the number of file outputs */
  number_of_outputs = 0;

  /* initialise the grid of locations */
  for (x = 0; x < grid_size; x++) {
    for (y = 0; y < grid_size; y++) {
      grid[x][y] = new Location(x, y);
    }
  }

  /* initialise the reaction diffusion IFN 1 system */
  rds.ifn_1 = new ReactDiff(RDType::ifn_1, status);
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* initialise the reaction diffusion IL 6 system */
  rds.il_6 = new ReactDiff(RDType::il_6, status);
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* initialise the reaction diffusion IL 10 system */
  rds.il_10 = new ReactDiff(RDType::il_10, status);
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* initialise the reaction diffusion IL 12 system */
  rds.il_12 = new ReactDiff(RDType::il_12, status);
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* initialise the reaction diffusion G-CSF system */
  rds.g_csf = new ReactDiff(RDType::g_csf, status);
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* initialise the reaction diffusion M-CSF system */
  rds.m_csf = new ReactDiff(RDType::m_csf, status);
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* initialise the reaction diffusion infected epithelial chemokine system */
  rds.epi_inf_ck = new ReactDiff(RDType::epi_inf_ck, status);
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* initialise the reaction diffusion NK apoptosis chemokine system */
  rds.apop_ck = new ReactDiff(RDType::apop_ck, status);
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* initialise the reaction diffusion Macrophage phagocytosis chemokine system */
  rds.phago_ck = new ReactDiff(RDType::phago_ck, status);
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* initialise reaction diffusion infectious virus model */
  rds.virus_inf = new ReactDiff(RDType::virus_inf, status);
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* initialise reaction diffusion non-infectious virus model */
  rds.virus_ninf = new ReactDiff(RDType::virus_ninf, status);
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* initialise the vessels in the tissue environment */
  (*status) = init_vessels();
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* initialise the epithelium */
  num_epi_in_initial_box = 0;
  (*status) = init_epithelium();
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* initialise the virions */
  (*status) = init_virions();
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* initialise the agents */
  (*status) = init_all_agents();
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  /* successful */
  (*status) = CAM_SUCCESS;
}


/*
 * =============================================================================
 * Environment destructor
 * =============================================================================
 */
Environment::~Environment ()
{
  int i = 0, j = 0;

  /* delete the agents */
  for (i = 0; i < (int) all_agents.size(); i++) {
    delete all_agents[i];
  }

  /* delete the epithelium */
  for (i = 0; i < (int) epithelium.size(); i++) {
    delete epithelium[i];
  }

  /* delete the vessels */
  for (i = 0; i < (int) vessels.size(); i++) {
    delete vessels[i];
  }

  /* delete the reaction-diffusion systems */
  if (rds.virus_ninf != nullptr) {
    delete rds.virus_ninf;
  }

  if (rds.virus_inf != nullptr) {
    delete rds.virus_inf;
  }

  if (rds.phago_ck != nullptr) {
    delete rds.phago_ck;
  }

  if (rds.apop_ck != nullptr) {
    delete rds.apop_ck;
  }

  if (rds.epi_inf_ck != nullptr) {
    delete rds.epi_inf_ck;
  }

  if (rds.m_csf != nullptr) {
    delete rds.m_csf;
  }

  if (rds.g_csf != nullptr) {
    delete rds.g_csf;
  }

  if (rds.il_12 != nullptr) {
    delete rds.il_12;
  }

  if (rds.il_10 != nullptr) {
    delete rds.il_10;
  }

  if (rds.il_6 != nullptr) {
    delete rds.il_6;
  }

  if (rds.ifn_1 != nullptr) {
    delete rds.ifn_1;
  }

  /* delete the grid locations */
  for (i = 0; i < grid_size; i++) {
    for (j = 0; j < grid_size; j++) {
      delete grid[i][j];
    }
  }
}


/*
 * =============================================================================
 * Initialise the vessels in the tissue environment
 * =============================================================================
 */
int
Environment::init_vessels ()
{
  int status = CAM_ERROR_MISC;

  /* initialise the blood vessels */
  if (params.env_fixed_blood_vessel_placement) {
    /* fixed blood vessels */

    for (std::pair<int,int> coord:params.fixed_blood_vessels) {
      if (grid[coord.first][coord.second]->is_vessel()) {
        CAM_ERROR("vessel already at location", CAM_ERROR_EXISTS);
      }

      status = add_vessel(coord.first, coord.second, VesselType::blood);
      CAM_ERROR_CHECK(status, CAM_SUCCESS);
    }
  } else {
    /* random distribution of blood vessels */

    int nves = 0;
    int attempts = 0;
    const int max_attempts = grid_size*grid_size;

    while (nves < params.env_initial_number_blood_vessels && attempts < max_attempts) {
      /* Random (x,y) coordinates */
      int x = distribution_0_grid_size(rng);
      int y = distribution_0_grid_size(rng);

      if (grid[x][y]->is_vessel()) {
        attempts++;
        continue;
      }

      status = add_vessel(x, y, VesselType::blood);
      CAM_ERROR_CHECK(status, CAM_SUCCESS);

      nves++;
      attempts = 0;
    }

    if (attempts >= max_attempts) {
      CAM_ERROR("unable to initialise the blood vessels", CAM_ERROR_VALUE);
    }
  }

  /* initialise the lymphatic vessels */
  if (params.env_fixed_lymph_vessel_placement) {
    /* fixed lymphatic vessels */

    for (std::pair<int,int> coord:params.fixed_lymph_vessels) {
      if (grid[coord.first][coord.second]->is_vessel()) {
        CAM_ERROR("vessel already at location", CAM_ERROR_EXISTS);
      }

      status = add_vessel(coord.first, coord.second, VesselType::lymph);
      CAM_ERROR_CHECK(status, CAM_SUCCESS);
    }
  } else {
    /* random distribution of lymph vessels */

    int nves = 0;
    int attempts = 0;
    const int max_attempts = grid_size*grid_size;

    while (nves < params.env_initial_number_lymph_vessels && attempts < max_attempts) {
      /* random (x,y) coords */
      int x = distribution_0_grid_size(rng);
      int y = distribution_0_grid_size(rng);

      if (grid[x][y]->is_vessel()) {
        attempts++;
        continue;
      }

      status = add_vessel(x, y, VesselType::lymph);
      CAM_ERROR_CHECK(status, CAM_SUCCESS);

      nves++;
      attempts = 0;
    }

    if (attempts >= max_attempts) {
      CAM_ERROR("unable to initialise lymph vessels", CAM_ERROR_VALUE);
    }
  }

  /* store the neighbours */
  for (int i = 0; i < (int) vessels.size(); i++) {
    std::pair<int,int> coords = vessels[i]->get_coordinates();

    /* neighbourhood depth 1 */
    std::vector<std::pair<int, int>> nbrs = cam_util_neighbourhood_moore(coords, 1);
    for (std::pair<int,int> nb:nbrs) {
      vessels[i]->immediate_moore.push_back(nb);
    }

    /* neighbourhood depth max */
    nbrs = cam_util_neighbourhood_moore(coords, params.env_neighbourhood_max_level);
    for (std::pair<int,int> nb:nbrs) {
      vessels[i]->neighbours.push_back(nb);
    }
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * add vessel to the tissue environment
 * =============================================================================
 */
int
Environment::add_vessel (const int x, const int y, const VesselType tp)
{
  /* initialise return status */
  int status = CAM_ERROR_MISC;

  /* add to extracellular vector */
  Vessel* ves = new Vessel(x, y, tp, &rds, &status);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  vessels.push_back(ves);

  /* set grid location as a particular vessel */
  grid[x][y]->set_vessel(ves);

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * initialise the epithelium
 * =============================================================================
 */
int
Environment::init_epithelium ()
{
  int x = 0, y = 0;
  int status = CAM_ERROR_MISC;

  const int midpoint = (grid_size - 1) / 2;
  const int Lbox = 50;
  for (x = 0; x < grid_size; x++) {
    for (y = 0; y < grid_size; y++) {
      if (grid[x][y]->is_vessel()) {
        continue;
      }

      const int max_coord = ( (abs(x - midpoint) > abs(y - midpoint)) ? abs(x - midpoint) : abs(y - midpoint) );
      if (max_coord <= Lbox)
        num_epi_in_initial_box++;

      status = add_epithelial(x, y, EpithelialType::epithelial, EpithelialState::healthy);
      CAM_ERROR_CHECK(status, CAM_SUCCESS);
    }
  }

  /* store the neighbours */
  for (int i = 0; i < (int) epithelium.size(); i++) {
    std::pair<int,int> coords = epithelium[i]->get_coordinates();

    /* neighbourhood depth 1 */
    std::vector<std::pair<int, int>> nbrs = cam_util_neighbourhood_moore(coords, 1);
    for (std::pair<int,int> nb:nbrs) {
      epithelium[i]->immediate_moore.push_back(nb);
    }

    /* neighbourhood depth max */
    nbrs = cam_util_neighbourhood_moore(coords, params.env_neighbourhood_max_level);
    for (std::pair<int,int> nb:nbrs) {
      epithelium[i]->neighbours.push_back(nb);
    }
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Add an epithelial cell to environment
 * =============================================================================
 */
int
Environment::add_epithelial (const int x, const int y, const EpithelialType tp, const EpithelialState st)
{
  /* initialise return status */
  int status = CAM_ERROR_MISC;

  /* add to epithelium */
  Epithelial* epi = new Epithelial(x, y, tp, st, &rds, &status);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  epithelium.push_back(epi);

  /* set the index in epithelium vector */
  epi->index = ((int) epithelium.size()) - 1;

  /* set grid location as an epithelial cell */
  grid[x][y]->set_epithelial(epi);

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Initialise the extracellular virions
 * =============================================================================
 */
int
Environment::init_virions ()
{
  double area = 0.0, ic = 0.0;

  if (params.env_fixed_virion_placement) {
    for (int i = 0; i < params.env_initial_number_virions; i++) {
      int j = ( (i < (int) params.fixed_virions.size()) ? i : (int) params.fixed_virions.size() - 1 );
      std::pair<int, int> coord = params.fixed_virions[j];

      if (grid[coord.first][coord.second]->is_vessel()) {
        CAM_ERROR("location is a vessel", CAM_ERROR_BADARGS);
      }

      area = grid[coord.first][coord.second]->get_epithelial()->get_max_area();
      ic = 1.0 / area;

      /* increase the virion density at location */
      rds.virus_inf->inc_conc(coord.first, coord.second, ic);
    }

  } else {
    double init_virs = 0.0;
    double nvir = 0;
    int attempts = 0;
    const int max_attempts = grid_size*grid_size;
    const int midpoint = (grid_size - 1) / 2;
    const int Lbox = 50;

    if (params.env_initial_moi < CAM_TOL) {
      init_virs = (double) params.env_initial_number_virions;
    } else {
      init_virs = params.env_initial_moi * ((double) num_epi_in_initial_box);
    }

    while ((nvir < init_virs) && (attempts < max_attempts)) {
      const int lo = ( ((midpoint - Lbox) < 0) ? 0 : (midpoint - Lbox) );
      const int hi = ( ((midpoint + Lbox) > (grid_size - 1)) ? grid_size - 1 : (midpoint + Lbox) );
      int x = cam_util_choose_random_index(lo, hi);
      int y = cam_util_choose_random_index(lo, hi);

      if (grid[x][y]->is_vessel()) {
        attempts++;
        continue;
      }

      area = grid[x][y]->get_epithelial()->get_max_area();
      ic = 1.0 / area;

      /* increase the virion density at location */
      rds.virus_inf->inc_conc(x, y, ic);

      /* Increment nvir and reset attempts, then continue */
      nvir += 1.0;
      attempts = 0;
    }

    if (attempts >= max_attempts) {
      CAM_ERROR("unable to place virions", CAM_ERROR_BADARGS);
    }
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Calculate local particle numbers from local density concentration
 * =============================================================================
 */
double
Environment::density_to_particle (Epithelial* epi, const bool infectious)
{
  /* coordinates */
  std::pair<int, int> xy = epi->get_coordinates();

  /* get the number density */
  double vir_inf_nd = rds.virus_inf->get_conc(xy.first, xy.second);
  double vir_ninf_nd = rds.virus_ninf->get_conc(xy.first, xy.second);

  double area = epi->get_max_area();

  /* calculate the number of virions at location */
  double vir_inf_num_d = vir_inf_nd * area;
  double vir_ninf_num_d = vir_ninf_nd * area;

  return ( (infectious) ? vir_inf_num_d : vir_ninf_num_d );
}


/*
 * =============================================================================
 * Initialise all of the agents
 * =============================================================================
 */
int
Environment::init_all_agents ()
{
  int status = CAM_ERROR_MISC;

  /* macrophages */
  status = init_agent<Macrophage, MacrophageType, MacrophageState>(
    params.env_fixed_macrophage_placement, params.fixed_macrophages, params.env_initial_number_macrophages,
    MacrophageType::alveolar, MacrophageState::null);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* nks */
  status = init_agent<NK, NKType, NKState>(
    params.env_fixed_nk_placement, params.fixed_nks, params.env_initial_number_nks,
    NKType::cytotoxic, NKState::null);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* neutrophils */
  status = init_agent<Neutrophil, NeutrophilType, NeutrophilState>(
    params.env_fixed_neutrophil_placement, params.fixed_neutrophils, params.env_initial_number_neutrophils,
    NeutrophilType::null, NeutrophilState::null);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* monocytes */
  status = init_agent<Monocyte, MonocyteType, MonocyteState>(
    params.env_fixed_monocyte_placement, params.fixed_monocytes, params.env_initial_number_monocytes,
    MonocyteType::null, MonocyteState::null);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Template routine to initialise the specific agents
 * =============================================================================
 */
template <class T, typename T_TYP, typename T_ST>
int
Environment::init_agent (const bool fixed, const std::vector<std::pair<int,int>> fixed_coords, const int initial_number, const T_TYP tp, const T_ST st)
{
  int status = CAM_ERROR_MISC;

  if (fixed) {
    for (std::pair<int, int> coord:fixed_coords) {
      if (grid[coord.first][coord.second]->is_vessel()) {
        CAM_ERROR("location is a vessel", CAM_ERROR_BADARGS);
      }

      status = add_agent<T, T_TYP, T_ST>(coord.first, coord.second, tp, st, false, true);
      CAM_ERROR_CHECK(status, CAM_SUCCESS);
    }

  } else {
    int nagents = 0;
    int max_attempts = grid_size*grid_size;
    int attempts = 0;

    while (nagents < initial_number && attempts < max_attempts) {
      int x = distribution_0_grid_size(rng);
      int y = distribution_0_grid_size(rng);

      if (grid[x][y]->is_vessel()) {
        attempts++;
        continue;
      }

      status = add_agent<T, T_TYP, T_ST>(x, y, tp, st, false, false);
      if (status < CAM_SUCCESS) {
        attempts++;
        continue;
      } else {
        CAM_ERROR_CHECK(status, CAM_SUCCESS);
      }

      /* Increment nvir and reset attempts, then continue */
      attempts = 0;
      nagents++;
    }

    if (attempts >= max_attempts) {
      CAM_ERROR("unable to place agents", CAM_ERROR_BADARGS);
    }
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Add a agent (templated type) to environment
 * =============================================================================
 */
template <class T, typename T_TYP, typename T_ST>
int
Environment::add_agent (const int x, const int y, const T_TYP tp, const T_ST st, const bool is_converted_agent, const bool use_output)
{
  /* initialise return status */
  int status = CAM_ERROR_MISC;

  /* get the epithelial */
  Epithelial* epi = grid[x][y]->get_epithelial();

  /* create agent */
  T* ag = new T(x, y, tp, st, &rds, use_output, &status);
  if (ag == nullptr) {
    CAM_ERROR("unable to create agent", CAM_ERROR_NOMEM);
  } else {
    CAM_ERROR_CHECK(status, CAM_SUCCESS);
  }

  /* check available area */
  double area_to_add = ag->get_area();
  bool successful = ( (is_converted_agent) ? true : epi->check_area(area_to_add) );
  if (!successful) {
    delete ag;

    return CAM_WARN_AREA;
  }

  /* add to vector */
  all_agents.push_back(ag);
  ag->index = ((int) all_agents.size()) - 1;

  /* add agent as resident */
  status = epi->add_resident(ag);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  ag->set_as_resident_on_cell(epi);

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Remove an agent from the environment
 * =============================================================================
 */
int
Environment::remove_agent (const int del_i)
{
  const int sz = (int) all_agents.size();
  Agent* ta = all_agents[del_i];

  /* coordinates and size */
  std::pair<int, int> xy = all_agents[del_i]->get_coordinates();

  /* remove resident */
  int status = grid[xy.first][xy.second]->get_epithelial()->remove_resident(ta);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* swap with end value of vector */
  all_agents[sz-1]->index = del_i;
  all_agents[del_i] = all_agents[sz-1];

  /* shrink vector */
  all_agents.pop_back();

  /* delete the agent (remove secretions etc..) */
  delete ta;

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Remove an agent from the environment
 * =============================================================================
 */
int
Environment::move_agent (const int agent_index, const int epithelial_index)
{
  int status = CAM_ERROR_MISC;

  /* coordinates */
  std::pair<int, int> xy = all_agents[agent_index]->get_coordinates();

  /* decrease agent cytokine secretion and uptakes */
  all_agents[agent_index]->dec_rds_values();

  /* remove resident */
  status = grid[xy.first][xy.second]->get_epithelial()->remove_resident(all_agents[agent_index]);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* add resident to new epithelial */
  status = epithelium[epithelial_index]->add_resident(all_agents[agent_index]);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* increase agent cytokine secretion and uptakes at new location */
  all_agents[agent_index]->inc_rds_values();

  /* re-define resident epithelial pointer in agent */
  all_agents[agent_index]->set_as_resident_on_cell(epithelium[epithelial_index]);

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Advance the time step
 * =============================================================================
 */
int
Environment::advance (const double time_n, const double time)
{
  int status = CAM_ERROR_MISC;

  /* update the age of all epithelial cells */
  for (int i = 0; i < (int) epithelium.size(); i++) {
    epithelium[i]->increase_age(timestep);
  }

  /* update the age of all agents */
  for (int i = 0; i < (int) all_agents.size(); i++) {
    all_agents[i]->increase_age(timestep);
  }


  /* compute the spread of the virus into the extracellular environment */
  status = compute_spread(time_n, time);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* compute the epithelial death */
  status = compute_epithelial_death(time);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* compute the agent actions */
  status = compute_agent_action(time);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* compute the agent migration */
  status = compute_agent_migration(time);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* compute the agent recruitment */
  status = compute_agent_recruitment();
  CAM_ERROR_CHECK(status, CAM_SUCCESS);


  /* advance the ifn_1 system */
  status = rds.ifn_1->advance();
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* advance the il_6 system */
  status = rds.il_6->advance();
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* advance the il_10 system */
  status = rds.il_10->advance();
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* advance the il_12 system */
  status = rds.il_12->advance();
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* advance the g_csf system */
  status = rds.g_csf->advance();
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* advance the m_csf system */
  status = rds.m_csf->advance();
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* advance the epi_inf_ck system */
  status = rds.epi_inf_ck->advance();
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* advance the NK apoptosis chemokine system */
  status = rds.apop_ck->advance();
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* advance the macrophage phagocytosis chemokine system */
  status = rds.phago_ck->advance();
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* advance the infectious virus system */
  status = rds.virus_inf->advance();
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  /* advance the non-infectious virus system */
  status = rds.virus_ninf->advance();
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Compute the viral spread into the extracellular environment
 * =============================================================================
 */
int
Environment::compute_spread (const double time_n, const double time)
{
  int status = CAM_ERROR_MISC;

  for (int i = 0; i < (int) epithelium.size(); i++) {
    /* update the extracellular virions vector */
    double num_virs = density_to_particle(epithelium[i], true);

    /* update receptor model and compute replication */
    status = epithelium[i]->replicate(num_virs, time_n, time);
    CAM_ERROR_CHECK(status, CAM_SUCCESS);
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Compute the death of the epithelial cells (bursting, apoptosis, etc..)
 * =============================================================================
 */
int
Environment::compute_epithelial_death (const double time)
{
  for (int i = 0; i < (int) epithelium.size(); i++) {
    Epithelial* epi = epithelium[i];

    /* coordinates */
    std::pair<int, int> xy = epi->get_coordinates();

    /* if not bursting, continue */
    if (epi->check_flag_bursting()) {
      std::pair<int, int> nbr;
      int ind = 0;

      /* scatter all epithelial residents */
      int j = 0;
      std::vector<Agent*>* residents = epi->get_residents();
      while (j < (int) (*residents).size()) { // size of residents is decreasing to empty
        /* using maximum area of epithelial cells, find a space */
        ind = find_space((*residents)[j], epi->neighbours);
        if (ind < 0) {
          /* no spaces available; remove agent */
          remove_agent((*residents)[j]->index);
          continue;
        }
        nbr = epi->neighbours[ind];

        /* move agent */
        move_agent((*residents)[j]->index, grid[nbr.first][nbr.second]->get_epithelial()->index);
      }

      /* generate random numbers */
      std::vector<double> rands;
      double r_sum = 0.0;
      for (int i = 0; i < (int) epi->neighbours.size(); i++) {
        double r = distribution_0_1(rng);
        rands.push_back(r);
        r_sum += r;
      }

      /* scatter epithelial interior */
      double conc = epi->get_interior() / epi->get_max_area();
      for (int i = 0; i < (int) epi->neighbours.size(); i++) {
        nbr = epi->neighbours[i];

        if (grid[nbr.first][nbr.second]->is_vessel()) {
          continue;
        }

        if (grid[nbr.first][nbr.second]->get_epithelial()->check_state(EpithelialState::burst)) {
          continue;
        }

        double r_conc = rands[i] * conc / r_sum;
        rds.virus_inf->inc_conc(nbr.first, nbr.second, r_conc);
      }

      /* re-generate random numbers */
      r_sum = 0.0;
      for (int i = 0; i < (int) epi->neighbours.size(); i++) {
        double r = distribution_0_1(rng);
        rands[i] = r;
        r_sum += r;
      }

      /* move the viral concentration */
      conc = rds.virus_inf->get_conc(xy.first, xy.second);
      for (int i = 0; i < (int) epi->neighbours.size(); i++) {
        nbr = epi->neighbours[i];

        if (grid[nbr.first][nbr.second]->is_vessel()) {
          continue;
        }

        if (grid[nbr.first][nbr.second]->get_epithelial()->check_state(EpithelialState::burst)) {
          continue;
        }

        double r_conc = rands[i] * conc / r_sum;
        rds.virus_inf->inc_conc(nbr.first, nbr.second, r_conc);
      }
      rds.virus_inf->set_decay(xy.first, xy.second, 1.0);
      rds.virus_inf->reset(xy.first, xy.second);

      /* set the epithelial cell as burst */
      epi->change_state(EpithelialState::burst);
      epi->switch_flag_bursting();
    }


    /* apoptosis */
    if (epi->check_flag_apoptosis()) {
      if (epi->is_being_apoptosed) {
        double pt = time - epi->apoptosis_start_time;
        if (pt > epi->apoptosis_timescale) {
          /* apoptosis is complete */
          epi->end_apoptosis();

          const int sz = (int) epi->neighbours.size();

          /* generate random numbers */
          std::vector<double> rands;
          double r_sum = 0.0;
          for (int i = 0; i < sz; i++) {
            double r = distribution_0_1(rng);
            rands.push_back(r);
            r_sum += r;
          }

          /* scatter epithelial interior */
          double conc = distribution_0_1(rng) * (epi->get_interior() / epi->get_max_area());
          for (int i = 0; i < sz; i++) {
            std::pair<int, int> nbr = epi->neighbours[i];

            if (grid[nbr.first][nbr.second]->is_vessel()) {
              continue;
            }

            if (grid[nbr.first][nbr.second]->get_epithelial()->check_state(EpithelialState::burst)) {
              continue;
            }

            double r_conc = rands[i] * conc / r_sum;
            rds.virus_inf->inc_conc(nbr.first, nbr.second, r_conc);
          }

          /* reset the virus rd components */
          rds.virus_inf->reset(xy.first, xy.second);
        }
      } else {
        /* start apoptosis from current time */
        epi->start_apoptosis(time);
      }
    }


    /* phagocytosis by macrophages & neutrophils */
    if (epi->check_flag_phagocytosis()) {
      if (epi->is_being_phagocytosed) {
        double pt = time - epi->phagocytosis_start_time;
        if (pt > epi->phagocytosis_timescale) {
          /* phagocytosis is complete */
          epi->end_phagocytosis();

          /* reset the virus rd components */
          rds.virus_inf->reset(xy.first, xy.second);
        }
      } else {
        /* start phagocytosis from current time */
        epi->start_phagocytosis(time);
      }
    }

  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Compute the death and actions of agent
 * =============================================================================
 */
int
Environment::compute_agent_action (const double time)
{
  int i = 0;
  while (i < (int) all_agents.size()) {
    if (all_agents[i]->exceeded_lifespan()) {
      /* if age of agent has exceeded its lifespan */

      if (all_agents[i]->is_being_apoptosed) {
        double pt = time - all_agents[i]->apoptosis_start_time;
        if (pt > all_agents[i]->apoptosis_timescale) {
          /* apoptosis is complete */
          all_agents[i]->end_apoptosis();

          if (all_agents[i]->is_phagocyte()) {
            Epithelial* epi = (Epithelial*) all_agents[i]->get_resident_on_cell();
            const int sz = (int) epi->neighbours.size();

            /* generate random numbers */
            std::vector<double> rands;
            double r_sum = 0.0;
            for (int i = 0; i < sz; i++) {
              double r = distribution_0_1(rng);
              rands.push_back(r);
              r_sum += r;
            }

            /* scatter epithelial interior */
            double conc = distribution_0_1(rng) * all_agents[i]->get_internalised();
            for (int i = 0; i < sz; i++) {
              std::pair<int, int> nbr = epi->neighbours[i];

              if (grid[nbr.first][nbr.second]->is_vessel()) {
                continue;
              }

              if (grid[nbr.first][nbr.second]->get_epithelial()->check_state(EpithelialState::burst)) {
                continue;
              }

              double r_conc = rands[i] * conc / r_sum;
              rds.virus_inf->inc_conc(nbr.first, nbr.second, r_conc);
            }
          }

        }
      } else {
        all_agents[i]->start_apoptosis(time);
      }

      i++;
      continue;
    }

    /* bursting */
    if (all_agents[i]->check_flag_bursting()) {
      std::pair<int, int> nbr;

      Epithelial* epi = (Epithelial*) all_agents[i]->get_resident_on_cell();

      /* generate random numbers */
      std::vector<double> rands;
      double r_sum = 0.0;
      for (int i = 0; i < (int) epi->neighbours.size(); i++) {
        double r = distribution_0_1(rng);
        rands.push_back(r);
        r_sum += r;
      }

      /* scatter agent interior */
      double conc = all_agents[i]->get_internalised();
      for (int i = 0; i < (int) epi->neighbours.size(); i++) {
        nbr = epi->neighbours[i];

        if (grid[nbr.first][nbr.second]->is_vessel()) {
          continue;
        }

        if (grid[nbr.first][nbr.second]->get_epithelial()->check_state(EpithelialState::burst)) {
          continue;
        }

        double r_conc = rands[i] * conc / r_sum;
        rds.virus_inf->inc_conc(nbr.first, nbr.second, r_conc);
      }

      /* remove the agent as it has burst */
      remove_agent(i);
      continue;
    }

    /* apoptosis */
    if (all_agents[i]->is_being_apoptosed) {
      double pt = time - all_agents[i]->apoptosis_start_time;
      if (pt > all_agents[i]->apoptosis_timescale) {
        /* apoptosis is complete */
        all_agents[i]->end_apoptosis();

        if (all_agents[i]->is_phagocyte()) {
          Epithelial* epi = (Epithelial*) all_agents[i]->get_resident_on_cell();
          const int sz = (int) epi->neighbours.size();

          /* generate random numbers */
          std::vector<double> rands;
          double r_sum = 0.0;
          for (int i = 0; i < sz; i++) {
            double r = distribution_0_1(rng);
            rands.push_back(r);
            r_sum += r;
          }

          /* scatter epithelial interior */
          double conc = distribution_0_1(rng) * all_agents[i]->get_internalised();
          for (int i = 0; i < sz; i++) {
            std::pair<int, int> nbr = epi->neighbours[i];

            if (grid[nbr.first][nbr.second]->is_vessel()) {
              continue;
            }

            if (grid[nbr.first][nbr.second]->get_epithelial()->check_state(EpithelialState::burst)) {
              continue;
            }

            double r_conc = rands[i] * conc / r_sum;
            rds.virus_inf->inc_conc(nbr.first, nbr.second, r_conc);
          }
        }

      }

      i++;
      continue;
    }

    /* phagocytosis */
    if (all_agents[i]->is_being_phagocytosed) {
      double pt = time - all_agents[i]->phagocytosis_start_time;
      if (pt > all_agents[i]->phagocytosis_timescale) {
        /* phagocytosis is complete */
        all_agents[i]->end_phagocytosis();

        remove_agent(i);
        continue;
      }

      i++;
      continue;
    }

    /* AGENT ACTION:
     *    macrophages (activation, deactivation and phagocytosis)
     *    neutrophils (activation, deactivation, virion consumption and phagocytosis)
     *    nks (activation, deactivation and trigger apoptosis)
     *    monocytes (activation, deactivation and conversion into macs)
     */
    all_agents[i]->action(time);

    /* check for monocyte conversion into macrophages */
    if (all_agents[i]->check_flag_conversion()) {
      /* converting monocyte into macrophage */

      /* coords */
      std::pair<int, int> xy = all_agents[i]->get_coordinates();

      /* add macrophage to extracellular environment */
      int status = add_agent<Macrophage, MacrophageType, MacrophageState>(
        xy.first, xy.second, MacrophageType::alveolar, MacrophageState::resting, true, false);
      CAM_ERROR_CHECK(status, CAM_SUCCESS);

      /* switch conversion flag */
      all_agents[i]->switch_flag_conversion();

      /* remove the monocyte; swaps with last entry in all_agents (new macrophage) */
      remove_agent(i);

      /* re-compute the action of the new macrophage to update secretions and uptakes */
      all_agents[i]->action(time);

      /* 'i' is incremented below as new agent is at position 'i' and therefore,
       * does not want to be computed this timestep */
    }

    /* increment i */
    i++;
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Compute the migration of the agents
 * =============================================================================
 */
int
Environment::compute_agent_migration (const double time)
{
  /* Migrate routines require the find_space and check_space routines from the
   * environment class. Therefore, we need to generate function pointers using
   * std::bind */
  auto fndsp = std::bind(&Environment::find_space, this, std::placeholders::_1, std::placeholders::_2);
  auto chksp = std::bind(&Environment::check_space, this, std::placeholders::_1, std::placeholders::_2);

  int i = 0;
  int attempts = 0;
  while (i < (int) all_agents.size()) {
    Epithelial* epi = (Epithelial*) all_agents[i]->get_resident_on_cell();

    /* find the index of the neighbour to move to */
    int ind = all_agents[i]->migrate(timestep, time, fndsp, chksp);
    if (ind < 0) {
      /* remove the agent, as migrated into vessel or empty space */
      remove_agent(i);
      continue;
    } else if (ind == ((int) epi->immediate_moore.size())) {
      /* not moving */
      i++;
      continue;
    }

    /* get the neighbour coords */
    std::pair<int, int> nbr = epi->immediate_moore[ind];

    if (grid[nbr.first][nbr.second]->get_epithelial()->is_being_phagocytosed) {
      attempts++;
      if (attempts > (int) epi->immediate_moore.size()) {
        /* not moving */
        i++;
        attempts = 0;
      }

      continue;
    }

    /* move the agent */
    move_agent(i, grid[nbr.first][nbr.second]->get_epithelial()->index);

    /* increment i */
    i++;
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Compute the recruitment of the agents
 * =============================================================================
 */
int
Environment::compute_agent_recruitment ()
{
  int status = CAM_ERROR_MISC;

  for (int i = 0; i < (int) vessels.size(); i++) {
    if (vessels[i]->get_type() != VesselType::blood) {
      /* can only recruit from blood vessels */
      continue;
    }

    /* recruit macrophages */
    status = recruit_agents<Macrophage, MacrophageType, MacrophageState>(vessels[i], MacrophageType::alveolar, MacrophageState::null);
    CAM_ERROR_CHECK(status, CAM_SUCCESS);

    /* recruit nks */
    status = recruit_agents<NK, NKType, NKState>(vessels[i], NKType::cytotoxic, NKState::null);
    CAM_ERROR_CHECK(status, CAM_SUCCESS);

    /* recruit neutrophils */
    status = recruit_agents<Neutrophil, NeutrophilType, NeutrophilState>(vessels[i], NeutrophilType::null, NeutrophilState::null);
    CAM_ERROR_CHECK(status, CAM_SUCCESS);

    /* recruit monocytes */
    status = recruit_agents<Monocyte, MonocyteType, MonocyteState>(vessels[i], MonocyteType::null, MonocyteState::null);
    CAM_ERROR_CHECK(status, CAM_SUCCESS);
  }

  return CAM_SUCCESS;
}


template <class T, typename T_TYP, typename T_ST>
int
Environment::recruit_agents (Vessel* ves, const T_TYP tp, const T_ST st)
{
  int status = CAM_ERROR_MISC;

  /* coordinates of vessel */
  std::pair<int, int> xy = ves->get_coordinates();

  /* create temporary agent so we can know the cytokine concs local to vessel
   * for recruitment of type T */
  T* temp = new T(xy.first, xy.second, tp, st, &rds, false, &status);
  if (temp == nullptr) {
    CAM_ERROR("unable to create temporary agent", CAM_ERROR_NOMEM);
  } else {
    CAM_ERROR_CHECK(status, CAM_SUCCESS);
  }

  /* inside the agent, check whether we want to recruit using specific cytokines */
  bool recruit_agent = temp->recruit(ves->immediate_moore);

  /* recruit */
  if (recruit_agent) {
    /* find a space */
    int ind = find_space(temp, ves->immediate_moore);
    if (!(ind < 0)) {
      std::pair<int, int> nbr = ves->immediate_moore[ind];

      /* add to extracellular environment */
      status = add_agent<T, T_TYP, T_ST>(nbr.first, nbr.second, tp, st, false, false);
      if (status > 0) {
        /* remove temporary memory */
        delete temp;

        CAM_ERROR("unble to recruit agent", CAM_ERROR_MISC);
      }
    }
  }

  /* remove temporary memory */
  delete temp;

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Find an available space in neighbours using maximum area of epithelial cell
 * =============================================================================
 */
int
Environment::find_space (Agent* agent_to_move, std::vector<std::pair<int,int>> neighbours)
{
  int index = -1;
  std::pair<int, int> xy;
  Epithelial* epi;

  std::vector<int> spaces;
  std::vector<Agent*> ags;
  for (int i = 0; i < (int) neighbours.size(); i++) {
    xy = neighbours[i];

    if (grid[xy.first][xy.second]->is_vessel()) {
      spaces.push_back(i);
      continue;
    }

    /* get the epithelial */
    epi = grid[xy.first][xy.second]->get_epithelial();
    if (epi->check_state(EpithelialState::burst) || epi->check_state(EpithelialState::phagocytosed)) {
      spaces.push_back(i);
      continue;
    }

    /* using the area, determine whether or not a move is possible */
    double area_to_add = agent_to_move->get_area();
    bool successful = epi->check_area(area_to_add);
    if (!successful) {
      /* if agent to move is inflammatory, check whether NK cells are resident
       * on space moving to, if so then we allow the move even though there isn't
       * space */
      if (agent_to_move->is_inflammatory() && epi->has_residents(AgentType::nk)) {
        spaces.push_back(i);
        continue;
      }

      /* if agent to move is NK cell and neighbours have inflam. immune cells,
       * then we allow the move even though there isn't space. */
      if ((agent_to_move->get_agent_type() == AgentType::nk) && epi->has_residents_inflammatory()) {
        spaces.push_back(i);
        continue;
      }

      /* if agent to move is phagocyte and neighbours have apoptotic immune cells
       * which are waiting to be phagocytosed then allow the move even though there
       * isn't space. */
      if (agent_to_move->is_phagocyte() && epi->has_residents_apoptotic()) {
        spaces.push_back(i);
        continue;
      }

      continue;
    }

    /* record the space */
    spaces.push_back(i);
  }

  if (!spaces.empty()) {
    /* The vector - spaces - contains all the spaces found in neighbours. We now
     * need to randomly choose a neighbour to return */
    int sp_i = cam_util_choose_random_index(0, ((int) spaces.size()) - 1);
    index = spaces[sp_i];

    /* check whether space is blood vessel */
    xy = neighbours[index];
    epi = grid[xy.first][xy.second]->get_epithelial();

    if (grid[xy.first][xy.second]->is_vessel()) {
      index = -2;
    } else if (epi->check_state(EpithelialState::burst)) {
      index = -3;
    } else if (epi->check_state(EpithelialState::phagocytosed)) {
      index = -4;
    }
  }

  return index;
}


/*
 * =============================================================================
 * Check whether the agent can move into a space
 * =============================================================================
 */
bool
Environment::check_space (Agent* agent_to_move, std::pair<int,int> xy)
{
  if (grid[xy.first][xy.second]->is_vessel()) {
    return false;
  }

  /* get the epithelial */
  Epithelial* epi = grid[xy.first][xy.second]->get_epithelial();
  if (epi == nullptr) {
    return false;
  } else if (epi->check_state(EpithelialState::burst) || epi->check_state(EpithelialState::phagocytosed) || epi->is_being_phagocytosed) {
    return false;
  }

  /* using the area, determine whether or not a move is possible */
  double area_to_add = agent_to_move->get_area();
  bool successful = epi->check_area(area_to_add);

  if (!successful) {
    /* if agent to move is inflammatory, check whether NK cells are resident
     * on space moving to, if so then we allow the move even though there isn't
     * space */
    if (agent_to_move->is_inflammatory() && epi->has_residents(AgentType::nk)) {
      return true;
    }

    /* if agent to move is NK cell and neighbours have inflam. immune cells,
     * then we allow the move even though there isn't space. */
    if ((agent_to_move->get_agent_type() == AgentType::nk) && epi->has_residents_inflammatory()) {
      return true;
    }

    /* if agent to move is phagocyte and neighbours have apoptotic immune cells
     * which are waiting to be phagocytosed then allow the move even though there
     * isn't space. */
    if (agent_to_move->is_phagocyte() && epi->has_residents_apoptotic()) {
      return true;
    }
  }

  return successful;
}


/*
 * =============================================================================
 * write the constants to the output file
 * =============================================================================
 */
void
Environment::write_consts (struct cam_output_files* files)
{
  /* seed */
  fprintf(files->consts, "%d\t", seed);

  /* number of outputs */
  fprintf(files->consts, "%d\t", number_of_outputs);

  /* simulation time */
  fprintf(files->consts, "%.3f\t", time_limit);

  /* time step */
  fprintf(files->consts, "%.3f\t", timestep);

  /* output interval */
  fprintf(files->consts, "%.3f\t", output_interval);

  /* grid size */
  fprintf(files->consts, "%d\n", grid_size);
}


/*
 * =============================================================================
 * write the data to the output files
 * =============================================================================
 */
void
Environment::write_data (struct cam_output_files* files, const double time)
{
  int x = 0, y = 0;


  /* initialise the grid codes */
  int** grid_codes = (int **) malloc(grid_size*sizeof(int*));
  for (x = 0; x < grid_size; x++) {
    grid_codes[x] = (int *) malloc(grid_size*sizeof(int));
  }
  CAM_OUTPUT_RESET_GRID_CODES(0);


  /* increment number of outputs */
  number_of_outputs++;


  /* For the first infected cell, plot intracellular replication data */
  if (params.env_fixed_virion_placement) {
    Epithelial* epi = grid[params.fixed_virions[0].first][params.fixed_virions[0].second]->get_epithelial();
    if (epi != nullptr) {
      epi->write_re_output(time);
      epi->write_vr_output(time);
    }
  }


  /*
   * ---------------------------------------------------------------------------
   * Location grid codes
   * ---------------------------------------------------------------------------
   */
  int num_epi_healthy = 0;
  int num_epi_infected = 0;
  int num_epi_apoptotic = 0;
  int num_epi_phagocytosed = 0;
  int num_epi_burst = 0;
  for (x = 0; x < grid_size; x++) {
    for (y = 0; y < grid_size; y++) {
      int grid_code = -1;

      if (grid[x][y]->is_vessel()) {
        grid_code = grid[x][y]->get_vessel()->get_code();
      } else if (grid[x][y]->is_epithelial()) {
        grid_code = grid[x][y]->get_epithelial()->get_code();

        num_epi_healthy += ( (grid[x][y]->get_epithelial()->check_state(EpithelialState::healthy)) ? 1 : 0 );
        num_epi_infected += ( (grid[x][y]->get_epithelial()->check_state(EpithelialState::infected)) ? 1 : 0 );
        num_epi_apoptotic += ( (grid[x][y]->get_epithelial()->check_state(EpithelialState::apoptotic)) ? 1 : 0 );
        num_epi_phagocytosed += ( (grid[x][y]->get_epithelial()->check_state(EpithelialState::phagocytosed)) ? 1 : 0 );
        num_epi_burst += ( (grid[x][y]->get_epithelial()->check_state(EpithelialState::burst)) ? 1 : 0 );
      }

      fprintf(files->grid_codes, "%d\t", grid_code);
    }
  }
  fprintf(files->grid_codes, "\n");

  fprintf(files->cell_counts_epithelial_healthy, "%.5e %d\n", time, num_epi_healthy);
  fprintf(files->cell_counts_epithelial_infected, "%.5e %d\n", time, num_epi_infected);
  fprintf(files->cell_counts_epithelial_apoptotic, "%.5e %d\n", time, num_epi_apoptotic);
  fprintf(files->cell_counts_epithelial_phagocytosed, "%.5e %d\n", time, num_epi_phagocytosed);
  fprintf(files->cell_counts_epithelial_burst, "%.5e %d\n", time, num_epi_burst);

  const int epi_tot = num_epi_infected + num_epi_apoptotic + num_epi_phagocytosed + num_epi_burst;
  const double epi_tot_d = ((double) epi_tot) / ((double) epithelium.size());
  fprintf(
    files->cell_counts_epithelial_ratio,
    "%.5e %g\n", time, epi_tot_d
  );


  /*
   * ---------------------------------------------------------------------------
   * Macrophage grid codes and counts
   * ---------------------------------------------------------------------------
   */

  /* resting alveolar macrophage */
  int num_macs_alv_rest = 0;
  CAM_OUTPUT_RESET_GRID_CODES(-1);
  num_macs_alv_rest = get_agent_codes<Macrophage, MacrophageType, MacrophageState>(
    grid_codes, MacrophageType::alveolar, MacrophageState::resting, time);
  CAM_OUTPUT_WRITE_GRID_CODES(files->grid_codes_macrophage_alveolar_resting);


  /* active alveolar macrophage */
  int num_macs_alv_act = 0;
  CAM_OUTPUT_RESET_GRID_CODES(-1);
  num_macs_alv_act = get_agent_codes<Macrophage, MacrophageType, MacrophageState>(
    grid_codes, MacrophageType::alveolar, MacrophageState::active, time);
  CAM_OUTPUT_WRITE_GRID_CODES(files->grid_codes_macrophage_alveolar_active);


  /* apoptotic alveolar macrophage */
  int num_macs_alv_apop = 0;
  CAM_OUTPUT_RESET_GRID_CODES(-1);
  num_macs_alv_apop = get_agent_codes<Macrophage, MacrophageType, MacrophageState>(
    grid_codes, MacrophageType::alveolar, MacrophageState::apoptotic, time);
  CAM_OUTPUT_WRITE_GRID_CODES(files->grid_codes_macrophage_alveolar_apoptotic);


  /* counts of macrophages */
  fprintf(files->cell_counts_macrophages_alveolar_resting, "%.5e %d\n", time, num_macs_alv_rest);
  fprintf(files->cell_counts_macrophages_alveolar_active, "%.5e %d\n", time, num_macs_alv_act);
  fprintf(files->cell_counts_macrophages_alveolar_apoptotic, "%.5e %d\n", time, num_macs_alv_apop);

  fprintf(files->cell_counts_macrophages_alveolar, "%.5e %d\n", time,
    num_macs_alv_rest + num_macs_alv_act + num_macs_alv_apop);


  /*
   * ---------------------------------------------------------------------------
   * NK grid codes and counts
   * ---------------------------------------------------------------------------
   */

  /* resting cytotoxic nks */
  int num_nks_cytotoxic_rest = 0;
  CAM_OUTPUT_RESET_GRID_CODES(-1);
  num_nks_cytotoxic_rest = get_agent_codes<NK, NKType, NKState>(
    grid_codes, NKType::cytotoxic, NKState::resting, time);
  CAM_OUTPUT_WRITE_GRID_CODES(files->grid_codes_nk_cytotoxic_resting);


  /* active cytotoxic nks */
  int num_nks_cytotoxic_act = 0;
  CAM_OUTPUT_RESET_GRID_CODES(-1);
  num_nks_cytotoxic_act = get_agent_codes<NK, NKType, NKState>(
    grid_codes, NKType::cytotoxic, NKState::active, time);
  CAM_OUTPUT_WRITE_GRID_CODES(files->grid_codes_nk_cytotoxic_active);


  /* apoptotic cytotoxic nks */
  int num_nks_cytotoxic_apop = 0;
  CAM_OUTPUT_RESET_GRID_CODES(-1);
  num_nks_cytotoxic_apop = get_agent_codes<NK, NKType, NKState>(
    grid_codes, NKType::cytotoxic, NKState::apoptotic, time);
  CAM_OUTPUT_WRITE_GRID_CODES(files->grid_codes_nk_cytotoxic_apoptotic);


  /* counts of nks */
  fprintf(files->cell_counts_nks_cytotoxic_resting, "%.5e %d\n", time, num_nks_cytotoxic_rest);
  fprintf(files->cell_counts_nks_cytotoxic_active, "%.5e %d\n", time, num_nks_cytotoxic_act);
  fprintf(files->cell_counts_nks_cytotoxic_apoptotic, "%.5e %d\n", time, num_nks_cytotoxic_apop);

  fprintf(files->cell_counts_nks_cytotoxic, "%.5e %d\n", time,
    num_nks_cytotoxic_rest + num_nks_cytotoxic_act + num_nks_cytotoxic_apop);


  /*
   * ---------------------------------------------------------------------------
   * Neutrophil grid codes and counts
   * ---------------------------------------------------------------------------
   */

  /* resting neutrophils */
  int num_neut_null_rest = 0;
  CAM_OUTPUT_RESET_GRID_CODES(-1);
  num_neut_null_rest = get_agent_codes<Neutrophil, NeutrophilType, NeutrophilState>(
    grid_codes, NeutrophilType::null, NeutrophilState::resting, time);
  CAM_OUTPUT_WRITE_GRID_CODES(files->grid_codes_neutrophil_null_resting);


  /* active neutrophils */
  int num_neut_null_act = 0;
  CAM_OUTPUT_RESET_GRID_CODES(-1);
  num_neut_null_act = get_agent_codes<Neutrophil, NeutrophilType, NeutrophilState>(
    grid_codes, NeutrophilType::null, NeutrophilState::active, time);
  CAM_OUTPUT_WRITE_GRID_CODES(files->grid_codes_neutrophil_null_active);

  /* apoptotic neutrophils */
  int num_neut_null_apop = 0;
  CAM_OUTPUT_RESET_GRID_CODES(-1);
  num_neut_null_apop = get_agent_codes<Neutrophil, NeutrophilType, NeutrophilState>(
    grid_codes, NeutrophilType::null, NeutrophilState::apoptotic, time);
  CAM_OUTPUT_WRITE_GRID_CODES(files->grid_codes_neutrophil_null_apoptotic);


  /* counts of neutrophils */
  fprintf(files->cell_counts_neutrophils_null_resting, "%.5e %d\n", time, num_neut_null_rest);
  fprintf(files->cell_counts_neutrophils_null_active, "%.5e %d\n", time, num_neut_null_act);
  fprintf(files->cell_counts_neutrophils_null_apoptotic, "%.5e %d\n", time, num_neut_null_apop);

  fprintf(files->cell_counts_neutrophils_null, "%.5e %d\n", time,
    num_neut_null_rest + num_neut_null_act + num_neut_null_apop);


  /*
   * ---------------------------------------------------------------------------
   * Monocyte grid codes and counts
   * ---------------------------------------------------------------------------
   */

  /* resting monocytes */
  int num_mono_null_rest = 0;
  CAM_OUTPUT_RESET_GRID_CODES(-1);
  num_mono_null_rest = get_agent_codes<Monocyte, MonocyteType, MonocyteState>(
    grid_codes, MonocyteType::null, MonocyteState::resting, time);
  CAM_OUTPUT_WRITE_GRID_CODES(files->grid_codes_monocyte_null_resting);


  /* active monocytes */
  int num_mono_null_act = 0;
  CAM_OUTPUT_RESET_GRID_CODES(-1);
  num_mono_null_act = get_agent_codes<Monocyte, MonocyteType, MonocyteState>(
    grid_codes, MonocyteType::null, MonocyteState::active, time);
  CAM_OUTPUT_WRITE_GRID_CODES(files->grid_codes_monocyte_null_active);


  /* apoptotic monocytes */
  int num_mono_null_apop = 0;
  CAM_OUTPUT_RESET_GRID_CODES(-1);
  num_mono_null_apop = get_agent_codes<Monocyte, MonocyteType, MonocyteState>(
    grid_codes, MonocyteType::null, MonocyteState::apoptotic, time);
  CAM_OUTPUT_WRITE_GRID_CODES(files->grid_codes_monocyte_null_apoptotic);


  /* counts of monocytes */
  fprintf(files->cell_counts_monocytes_null_resting, "%.5e %d\n", time, num_mono_null_rest);
  fprintf(files->cell_counts_monocytes_null_active, "%.5e %d\n", time, num_mono_null_act);
  fprintf(files->cell_counts_monocytes_null_apoptotic, "%.5e %d\n", time, num_mono_null_apop);

  fprintf(files->cell_counts_monocytes_null, "%.5e %d\n", time,
    num_mono_null_rest + num_mono_null_act + num_mono_null_apop);


  /*
   * ---------------------------------------------------------------------------
   * Reaction diffusion output and total cytokine
   * ---------------------------------------------------------------------------
   */
  double tot_ifn_1 = 0.0;
  double tot_il_6 = 0.0, tot_il_10 = 0.0, tot_il_12 = 0.0;
  double tot_g_csf = 0.0, tot_m_csf = 0.0;
  double tot_epi_inf_ck = 0.0, tot_apop_ck = 0.0, tot_phago_ck = 0.0;
  double tot_virus_inf = 0.0, tot_virus_ninf = 0.0;
  for (int x = 0; x < grid_size; x++) {
    for (int y = 0; y < grid_size; y++) {
      double tmp_area = 0.0;
      if (grid[x][y]->is_vessel()) {
        tmp_area = grid[x][y]->get_vessel()->get_max_area();
      } else if (grid[x][y]->is_epithelial()) {
        tmp_area = grid[x][y]->get_epithelial()->get_max_area();
      }

      /* cytokine values at location */
      double rdifn1 = ( (rds.ifn_1 == nullptr) ? 0.0 : rds.ifn_1->get_conc(x, y) );
      double rdil6 = ( (rds.il_6 == nullptr) ? 0.0 : rds.il_6->get_conc(x, y) );
      double rdil10 = ( (rds.il_10 == nullptr) ? 0.0 : rds.il_10->get_conc(x, y) );
      double rdil12 = ( (rds.il_12 == nullptr) ? 0.0 : rds.il_12->get_conc(x, y) );
      double rdgcsf = ( (rds.g_csf == nullptr) ? 0.0 : rds.g_csf->get_conc(x, y) );
      double rdmcsf = ( (rds.m_csf == nullptr) ? 0.0 : rds.m_csf->get_conc(x, y) );
      double rdepiinfck = ( (rds.epi_inf_ck == nullptr) ? 0.0 : rds.epi_inf_ck->get_conc(x, y) );
      double rdapopck = ( (rds.apop_ck == nullptr) ? 0.0 : rds.apop_ck->get_conc(x, y) );
      double rdphagock = ( (rds.phago_ck == nullptr) ? 0.0 : rds.phago_ck->get_conc(x, y) );
      double rdvrinf = ( (rds.virus_inf == nullptr) ? 0.0 : rds.virus_inf->get_conc(x, y) );
      double rdvrninf = ( (rds.virus_ninf == nullptr) ? 0.0 : rds.virus_ninf->get_conc(x, y) );

      /* total cytokine */
      tot_ifn_1 += rdifn1 * tmp_area;
      tot_il_6 += rdil6 * tmp_area;
      tot_il_10 += rdil10 * tmp_area;
      tot_il_12 += rdil12 * tmp_area;
      tot_g_csf += rdgcsf * tmp_area;
      tot_m_csf += rdmcsf * tmp_area;
      tot_epi_inf_ck += rdepiinfck * tmp_area;
      tot_apop_ck += rdapopck * tmp_area;
      tot_phago_ck += rdphagock * tmp_area;
      tot_virus_inf += rdvrinf * tmp_area;
      tot_virus_ninf += rdvrninf * tmp_area;

      fprintf(files->rd_ifn_1, "%.5e\t", rdifn1);
      fprintf(files->rd_il_6, "%.5e\t", rdil6);
      fprintf(files->rd_il_10, "%.5e\t", rdil10);
      fprintf(files->rd_il_12, "%.5e\t", rdil12);
      fprintf(files->rd_g_csf, "%.5e\t", rdgcsf);
      fprintf(files->rd_m_csf, "%.5e\t", rdmcsf);
      fprintf(files->rd_epi_inf_ck, "%.5e\t", rdepiinfck);
      fprintf(files->rd_apop_ck, "%.5e\t", rdapopck);
      fprintf(files->rd_phago_ck, "%.5e\t", rdphagock);
      fprintf(files->rd_virus_inf, "%.5e\t", rdvrinf);
      fprintf(files->rd_virus_ninf, "%.5e\t", rdvrninf);
    }
  }
  fprintf(files->rd_ifn_1, "\n");
  fprintf(files->rd_il_6, "\n");
  fprintf(files->rd_il_10, "\n");
  fprintf(files->rd_il_12, "\n");
  fprintf(files->rd_g_csf, "\n");
  fprintf(files->rd_m_csf, "\n");
  fprintf(files->rd_epi_inf_ck, "\n");
  fprintf(files->rd_apop_ck, "\n");
  fprintf(files->rd_phago_ck, "\n");
  fprintf(files->rd_virus_inf, "\n");
  fprintf(files->rd_virus_ninf, "\n");

  /* total cytokine */
  fprintf(files->rd_total_ifn_1, "%.5e %g\n", time, tot_ifn_1);
  fprintf(files->rd_total_il_6, "%.5e %g\n", time, tot_il_6);
  fprintf(files->rd_total_il_10, "%.5e %g\n", time, tot_il_10);
  fprintf(files->rd_total_il_12, "%.5e %g\n", time, tot_il_12);
  fprintf(files->rd_total_g_csf, "%.5e %g\n", time, tot_g_csf);
  fprintf(files->rd_total_m_csf, "%.5e %g\n", time, tot_m_csf);
  fprintf(files->rd_total_epi_inf_ck, "%.5e %g\n", time, tot_epi_inf_ck);
  fprintf(files->rd_total_apop_ck, "%.5e %g\n", time, tot_apop_ck);
  fprintf(files->rd_total_phago_ck, "%.5e %g\n", time, tot_phago_ck);
  fprintf(files->rd_total_virus_inf, "%.5e %g\n", time, tot_virus_inf);
  fprintf(files->rd_total_virus_ninf, "%.5e %g\n", time, tot_virus_ninf);
  fprintf(files->rd_total_virus, "%.5e %g\n", time, tot_virus_ninf + tot_virus_inf);


  /*
   * ---------------------------------------------------------------------------
   * Virus counts
   * ---------------------------------------------------------------------------
   */
  double num_virs_inf_extracellular = 0.0, num_virs_inf_intracellular = 0.0;
  double num_virs_ninf_extracellular = 0.0;
  for (int i = 0; i < (int) epithelium.size(); i++) {
    num_virs_inf_extracellular += density_to_particle(epithelium[i], true);
    num_virs_inf_intracellular += epithelium[i]->get_interior();

    num_virs_ninf_extracellular += density_to_particle(epithelium[i], false);
  }
  fprintf(files->cell_counts_virus_inf_extracellular, "%.5e %g\n", time, num_virs_inf_extracellular);
  fprintf(files->cell_counts_virus_inf_intracellular, "%.5e %g\n", time, num_virs_inf_intracellular);
  fprintf(files->cell_counts_virus_inf_total, "%.5e %g\n", time, num_virs_inf_extracellular + num_virs_inf_intracellular);

  fprintf(files->cell_counts_virus_ninf_extracellular, "%.5e %g\n", time, num_virs_ninf_extracellular);
  fprintf(files->cell_counts_virus_ninf_total, "%.5e %g\n", time, num_virs_ninf_extracellular);

  fprintf(files->cell_counts_virus_extracellular_total, "%.5e %g\n", time, num_virs_inf_extracellular + num_virs_ninf_extracellular);
  fprintf(files->cell_counts_virus_total, "%.5e %g\n", time, num_virs_inf_extracellular + num_virs_ninf_extracellular + num_virs_inf_intracellular);


  /* delete the grid_codes */
  for (int x = 0; x < grid_size; x++) {
    free(grid_codes[x]);
    grid_codes[x] = nullptr;
  }
  free(grid_codes);
  grid_codes = nullptr;
}


/*
 * =============================================================================
 * template routine to get the agent output codes
 * =============================================================================
 */
template <class T, typename T_TYP, typename T_ST>
int
Environment::get_agent_codes (int** grid_codes, const T_TYP tp, const T_ST st, const double time)
{
  int count = 0;

  for (int i = 0; i < (int) all_agents.size(); i++) {
    T* ag = dynamic_cast<T*>(all_agents[i]);
    if (ag == nullptr) {
      /* failed cast */
      continue;
    }

    /* coordinates */
    std::pair<int,int> xy = ag->get_coordinates();

    if ((ag->get_type() == tp) && (ag->get_state() == st)) {
      grid_codes[xy.first][xy.second] = ag->get_code();
      count++;

      ag->write_output(time);
    }
  }

  return count;
}
