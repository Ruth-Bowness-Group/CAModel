#include <cstring>

#include "cam_macrophage.h"
#include "macrophage.h"
#include "../../cam_error.h"
#include "../../cam_random.h"
#include "../../cam_util.h"
#include "../../cam_parameters.h"


/* struct containing macrophage specific parameters */
static macrophage_params params;


/*
 * =============================================================================
 * Set the macrophage specific parameters
 * =============================================================================
 */
int
cam_macrophage_set_params (YAML::Node* cfg)
{
  std::string ren;

  try {
    /* radius */
    params.macrophage_radius = (*cfg)["macrophage_radius"].as<double>();

    /* probability of activating macrophage */
    params.macrophage_active_probability = (*cfg)["macrophage_active_probability"].as<double>();

    /* lifespan */
    params.macrophage_resting_lifespan = (*cfg)["macrophage_resting_lifespan"].as<double>();
    params.macrophage_active_lifespan = (*cfg)["macrophage_active_lifespan"].as<double>();

    params.macrophage_lifespan_downreg_virus_EC50 = (*cfg)["macrophage_lifespan_downreg_virus_EC50"].as<double>();
    params.macrophage_lifespan_downreg_virus_hill = (*cfg)["macrophage_lifespan_downreg_virus_hill"].as<double>();

    /*
     * cytokines
     */

    /* T1 IFN */
    params.macrophage_resting_basal_secreted_ifn_1 = (*cfg)["macrophage_resting_basal_secreted_ifn_1"].as<double>();
    params.macrophage_resting_secreted_ifn_1 = (*cfg)["macrophage_resting_secreted_ifn_1"].as<double>();
    params.macrophage_resting_uptake_ifn_1 = (*cfg)["macrophage_resting_uptake_ifn_1"].as<double>();

    params.macrophage_active_basal_secreted_ifn_1 = (*cfg)["macrophage_active_basal_secreted_ifn_1"].as<double>();
    params.macrophage_active_secreted_ifn_1 = (*cfg)["macrophage_active_secreted_ifn_1"].as<double>();
    params.macrophage_active_uptake_ifn_1 = (*cfg)["macrophage_active_uptake_ifn_1"].as<double>();

    params.macrophage_ifn_1_EC50 = (*cfg)["macrophage_ifn_1_EC50"].as<double>();
    params.macrophage_ifn_1_hill = (*cfg)["macrophage_ifn_1_hill"].as<double>();
    params.macrophage_ifn_1_signal_virus_EC50 = (*cfg)["macrophage_ifn_1_signal_virus_EC50"].as<double>();
    params.macrophage_ifn_1_signal_virus_hill = (*cfg)["macrophage_ifn_1_signal_virus_hill"].as<double>();
    params.macrophage_ifn_1_downreg_il_10_EC50 = (*cfg)["macrophage_ifn_1_downreg_il_10_EC50"].as<double>();
    params.macrophage_ifn_1_downreg_il_10_hill = (*cfg)["macrophage_ifn_1_downreg_il_10_hill"].as<double>();

    /* IL-6 */
    params.macrophage_resting_basal_secreted_il_6 = (*cfg)["macrophage_resting_basal_secreted_il_6"].as<double>();
    params.macrophage_resting_secreted_il_6 = (*cfg)["macrophage_resting_secreted_il_6"].as<double>();
    params.macrophage_resting_uptake_il_6 = (*cfg)["macrophage_resting_uptake_il_6"].as<double>();

    params.macrophage_active_basal_secreted_il_6 = (*cfg)["macrophage_active_basal_secreted_il_6"].as<double>();
    params.macrophage_active_secreted_il_6 = (*cfg)["macrophage_active_secreted_il_6"].as<double>();
    params.macrophage_active_uptake_il_6 = (*cfg)["macrophage_active_uptake_il_6"].as<double>();

    params.macrophage_il_6_EC50 = (*cfg)["macrophage_il_6_EC50"].as<double>();
    params.macrophage_il_6_hill = (*cfg)["macrophage_il_6_hill"].as<double>();
    params.macrophage_il_6_signal_virus_EC50 = (*cfg)["macrophage_il_6_signal_virus_EC50"].as<double>();
    params.macrophage_il_6_signal_virus_hill = (*cfg)["macrophage_il_6_signal_virus_hill"].as<double>();
    params.macrophage_il_6_downreg_il_10_EC50 = (*cfg)["macrophage_il_6_downreg_il_10_EC50"].as<double>();
    params.macrophage_il_6_downreg_il_10_hill = (*cfg)["macrophage_il_6_downreg_il_10_hill"].as<double>();

    /* IL-10 */
    params.macrophage_resting_basal_secreted_il_10 = (*cfg)["macrophage_resting_basal_secreted_il_10"].as<double>();
    params.macrophage_resting_secreted_il_10 = (*cfg)["macrophage_resting_secreted_il_10"].as<double>();
    params.macrophage_resting_uptake_il_10 = (*cfg)["macrophage_resting_uptake_il_10"].as<double>();

    params.macrophage_active_basal_secreted_il_10 = (*cfg)["macrophage_active_basal_secreted_il_10"].as<double>();
    params.macrophage_active_secreted_il_10 = (*cfg)["macrophage_active_secreted_il_10"].as<double>();
    params.macrophage_active_uptake_il_10 = (*cfg)["macrophage_active_uptake_il_10"].as<double>();

    params.macrophage_il_10_EC50 = (*cfg)["macrophage_il_10_EC50"].as<double>();
    params.macrophage_il_10_hill = (*cfg)["macrophage_il_10_hill"].as<double>();
    params.macrophage_il_10_signal_il_6_EC50 = (*cfg)["macrophage_il_10_signal_il_6_EC50"].as<double>();
    params.macrophage_il_10_signal_il_6_hill = (*cfg)["macrophage_il_10_signal_il_6_hill"].as<double>();
    params.macrophage_il_10_signal_il_12_EC50 = (*cfg)["macrophage_il_10_signal_il_12_EC50"].as<double>();
    params.macrophage_il_10_signal_il_12_hill = (*cfg)["macrophage_il_10_signal_il_12_hill"].as<double>();
    params.macrophage_il_10_signal_g_csf_EC50 = (*cfg)["macrophage_il_10_signal_g_csf_EC50"].as<double>();
    params.macrophage_il_10_signal_g_csf_hill = (*cfg)["macrophage_il_10_signal_g_csf_hill"].as<double>();
    params.macrophage_il_10_signal_m_csf_EC50 = (*cfg)["macrophage_il_10_signal_m_csf_EC50"].as<double>();
    params.macrophage_il_10_signal_m_csf_hill = (*cfg)["macrophage_il_10_signal_m_csf_hill"].as<double>();

    /* IL-12 */
    params.macrophage_resting_basal_secreted_il_12 = (*cfg)["macrophage_resting_basal_secreted_il_12"].as<double>();
    params.macrophage_resting_secreted_il_12 = (*cfg)["macrophage_resting_secreted_il_12"].as<double>();
    params.macrophage_resting_uptake_il_12 = (*cfg)["macrophage_resting_uptake_il_12"].as<double>();

    params.macrophage_active_basal_secreted_il_12 = (*cfg)["macrophage_active_basal_secreted_il_12"].as<double>();
    params.macrophage_active_secreted_il_12 = (*cfg)["macrophage_active_secreted_il_12"].as<double>();
    params.macrophage_active_uptake_il_12 = (*cfg)["macrophage_active_uptake_il_12"].as<double>();

    params.macrophage_il_12_EC50 = (*cfg)["macrophage_il_12_EC50"].as<double>();
    params.macrophage_il_12_hill = (*cfg)["macrophage_il_12_hill"].as<double>();
    params.macrophage_il_12_signal_virus_EC50 = (*cfg)["macrophage_il_12_signal_virus_EC50"].as<double>();
    params.macrophage_il_12_signal_virus_hill = (*cfg)["macrophage_il_12_signal_virus_hill"].as<double>();
    params.macrophage_il_12_downreg_il_10_EC50 = (*cfg)["macrophage_il_12_downreg_il_10_EC50"].as<double>();
    params.macrophage_il_12_downreg_il_10_hill = (*cfg)["macrophage_il_12_downreg_il_10_hill"].as<double>();

    /* G-CSF */
    params.macrophage_resting_basal_secreted_g_csf = (*cfg)["macrophage_resting_basal_secreted_g_csf"].as<double>();
    params.macrophage_resting_secreted_g_csf = (*cfg)["macrophage_resting_secreted_g_csf"].as<double>();
    params.macrophage_resting_uptake_g_csf = (*cfg)["macrophage_resting_uptake_g_csf"].as<double>();

    params.macrophage_active_basal_secreted_g_csf = (*cfg)["macrophage_active_basal_secreted_g_csf"].as<double>();
    params.macrophage_active_secreted_g_csf = (*cfg)["macrophage_active_secreted_g_csf"].as<double>();
    params.macrophage_active_uptake_g_csf = (*cfg)["macrophage_active_uptake_g_csf"].as<double>();

    params.macrophage_g_csf_EC50 = (*cfg)["macrophage_g_csf_EC50"].as<double>();
    params.macrophage_g_csf_hill = (*cfg)["macrophage_g_csf_hill"].as<double>();
    params.macrophage_g_csf_signal_virus_EC50 = (*cfg)["macrophage_g_csf_signal_virus_EC50"].as<double>();
    params.macrophage_g_csf_signal_virus_hill = (*cfg)["macrophage_g_csf_signal_virus_hill"].as<double>();
    params.macrophage_g_csf_downreg_il_10_EC50 = (*cfg)["macrophage_g_csf_downreg_il_10_EC50"].as<double>();
    params.macrophage_g_csf_downreg_il_10_hill = (*cfg)["macrophage_g_csf_downreg_il_10_hill"].as<double>();

    /* M-CSF */
    params.macrophage_resting_basal_secreted_m_csf = (*cfg)["macrophage_resting_basal_secreted_m_csf"].as<double>();
    params.macrophage_resting_secreted_m_csf = (*cfg)["macrophage_resting_secreted_m_csf"].as<double>();
    params.macrophage_resting_uptake_m_csf = (*cfg)["macrophage_resting_uptake_m_csf"].as<double>();

    params.macrophage_active_basal_secreted_m_csf = (*cfg)["macrophage_active_basal_secreted_m_csf"].as<double>();
    params.macrophage_active_secreted_m_csf = (*cfg)["macrophage_active_secreted_m_csf"].as<double>();
    params.macrophage_active_uptake_m_csf = (*cfg)["macrophage_active_uptake_m_csf"].as<double>();

    params.macrophage_m_csf_EC50 = (*cfg)["macrophage_m_csf_EC50"].as<double>();
    params.macrophage_m_csf_hill = (*cfg)["macrophage_m_csf_hill"].as<double>();
    params.macrophage_m_csf_signal_virus_EC50 = (*cfg)["macrophage_m_csf_signal_virus_EC50"].as<double>();
    params.macrophage_m_csf_signal_virus_hill = (*cfg)["macrophage_m_csf_signal_virus_hill"].as<double>();
    params.macrophage_m_csf_downreg_il_10_EC50 = (*cfg)["macrophage_m_csf_downreg_il_10_EC50"].as<double>();
    params.macrophage_m_csf_downreg_il_10_hill = (*cfg)["macrophage_m_csf_downreg_il_10_hill"].as<double>();

    /* infected epithelial chemokine */
    params.macrophage_resting_uptake_epi_inf_ck = (*cfg)["macrophage_resting_uptake_epi_inf_ck"].as<double>();
    params.macrophage_active_uptake_epi_inf_ck = (*cfg)["macrophage_active_uptake_epi_inf_ck"].as<double>();

    /* apoptosis chemokine */
    params.macrophage_secreted_apop_ck = (*cfg)["macrophage_secreted_apop_ck"].as<double>();

    /* phagocytosis chemokine */
    params.macrophage_resting_uptake_phago_ck = (*cfg)["macrophage_resting_uptake_phago_ck"].as<double>();
    params.macrophage_active_uptake_phago_ck = (*cfg)["macrophage_active_uptake_phago_ck"].as<double>();

    params.macrophage_apoptotic_secreted_phago_ck = (*cfg)["macrophage_apoptotic_secreted_phago_ck"].as<double>();

    /* infectious virus */
    params.macrophage_resting_uptake_virus_inf = (*cfg)["macrophage_resting_uptake_virus_inf"].as<double>();
    params.macrophage_resting_decay_virus_inf = (*cfg)["macrophage_resting_decay_virus_inf"].as<double>();

    params.macrophage_active_uptake_virus_inf = (*cfg)["macrophage_active_uptake_virus_inf"].as<double>();
    params.macrophage_active_decay_virus_inf = (*cfg)["macrophage_active_decay_virus_inf"].as<double>();

    params.macrophage_virus_inf_upreg_extracellular_virus_EC50 = (*cfg)["macrophage_virus_inf_upreg_extracellular_virus_EC50"].as<double>();
    params.macrophage_virus_inf_upreg_extracellular_virus_hill = (*cfg)["macrophage_virus_inf_upreg_extracellular_virus_hill"].as<double>();
    params.macrophage_virus_inf_downreg_internalised_virus_EC50 = (*cfg)["macrophage_virus_inf_downreg_internalised_virus_EC50"].as<double>();
    params.macrophage_virus_inf_downreg_internalised_virus_hill = (*cfg)["macrophage_virus_inf_downreg_internalised_virus_hill"].as<double>();

    /* non-infectious virus - no global memory; not used */

    /* movement */
    params.macrophage_resting_movement_rate = (*cfg)["macrophage_resting_movement_rate"].as<double>();
    params.macrophage_resting_movement_bias = (*cfg)["macrophage_resting_movement_bias"].as<double>();

    params.macrophage_active_movement_rate = (*cfg)["macrophage_active_movement_rate"].as<double>();
    params.macrophage_active_movement_bias = (*cfg)["macrophage_active_movement_bias"].as<double>();

    params.macrophage_movement_upreg_virus_EC50 = (*cfg)["macrophage_movement_upreg_virus_EC50"].as<double>();
    params.macrophage_movement_upreg_virus_hill = (*cfg)["macrophage_movement_upreg_virus_hill"].as<double>();

    params.macrophage_movement_weight_epi_inf_ck = (*cfg)["macrophage_movement_weight_epi_inf_ck"].as<double>();
    params.macrophage_movement_weight_phago_ck = (*cfg)["macrophage_movement_weight_phago_ck"].as<double>();

    /* recruitment */
    params.macrophage_recruitment_weight_il_6 = (*cfg)["macrophage_recruitment_weight_il_6"].as<double>();
    params.macrophage_recruitment_weight_il_10 = (*cfg)["macrophage_recruitment_weight_il_10"].as<double>();
    params.macrophage_recruitment_weight_m_csf = (*cfg)["macrophage_recruitment_weight_m_csf"].as<double>();
    params.macrophage_recruitment_weight_epi_inf_ck = (*cfg)["macrophage_recruitment_weight_epi_inf_ck"].as<double>();
    params.macrophage_recruitment_weight_phago_ck = (*cfg)["macrophage_recruitment_weight_phago_ck"].as<double>();

    params.macrophage_recruitment_EC50 = (*cfg)["macrophage_recruitment_EC50"].as<double>();
    params.macrophage_recruitment_hill = (*cfg)["macrophage_recruitment_hill"].as<double>();

    /* activation and deactivation */
    params.macrophage_activate_weight_ifn_1 = (*cfg)["macrophage_activate_weight_ifn_1"].as<double>();
    params.macrophage_activate_weight_il_6 = (*cfg)["macrophage_activate_weight_il_6"].as<double>();
    params.macrophage_activate_weight_il_10 = (*cfg)["macrophage_activate_weight_il_10"].as<double>();
    params.macrophage_activate_weight_m_csf = (*cfg)["macrophage_activate_weight_m_csf"].as<double>();
    params.macrophage_activate_weight_epi_inf_ck = (*cfg)["macrophage_activate_weight_epi_inf_ck"].as<double>();
    params.macrophage_activate_weight_phago_ck = (*cfg)["macrophage_activate_weight_phago_ck"].as<double>();
    params.macrophage_activate_weight_virus_inf = (*cfg)["macrophage_activate_weight_virus_inf"].as<double>();

    params.macrophage_activate_EC50 = (*cfg)["macrophage_activate_EC50"].as<double>();
    params.macrophage_activate_hill = (*cfg)["macrophage_activate_hill"].as<double>();

    /* bursting */
    params.macrophage_resting_burst_threshold = (*cfg)["macrophage_resting_burst_threshold"].as<double>();
    params.macrophage_active_burst_threshold = (*cfg)["macrophage_active_burst_threshold"].as<double>();

    /* apoptosis */
    params.macrophage_apoptosis_timescale = (*cfg)["macrophage_apoptosis_timescale"].as<double>();

    /*==========*/
    /* not used */
    params.macrophage_apoptosis_max_rate = (*cfg)["macrophage_apoptosis_max_rate"].as<double>();
    params.macrophage_apoptosis_half_max = (*cfg)["macrophage_apoptosis_half_max"].as<double>();
    params.macrophage_apoptosis_hill = (*cfg)["macrophage_apoptosis_hill"].as<double>();
    /*==========*/

    params.macrophage_apoptosis_by_ifn_max_rate = (*cfg)["macrophage_apoptosis_by_ifn_max_rate"].as<double>();
    params.macrophage_apoptosis_by_ifn_half_max = (*cfg)["macrophage_apoptosis_by_ifn_half_max"].as<double>();
    params.macrophage_apoptosis_by_ifn_hill = (*cfg)["macrophage_apoptosis_by_ifn_hill"].as<double>();

    /* phagocytosis */
    params.macrophage_phagocytosis_timescale = (*cfg)["macrophage_phagocytosis_timescale"].as<double>();
    params.macrophage_active_phagocytose_apoptotic_epithelial_probability = (*cfg)["macrophage_active_phagocytose_apoptotic_epithelial_probability"].as<double>();
    params.macrophage_active_phagocytose_infected_epithelial_probability = (*cfg)["macrophage_active_phagocytose_infected_epithelial_probability"].as<double>();
    params.macrophage_active_phagocytose_immune_probability = (*cfg)["macrophage_active_phagocytose_immune_probability"].as<double>();

    params.macrophage_blocking_function = false;

  } catch (std::exception &e) {
    CAM_ERROR(e.what(), CAM_ERROR_IO);
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Macrophage constructor
 * =============================================================================
 */
Macrophage::Macrophage (const int x, const int y, const MacrophageType tp, const MacrophageState st,
                        struct cam_rd_sys* rds, const bool use_output, int* status)
  : Agent{ x, y, AgentType::macrophage, rds }
{
  /* agent is constructed first */

  /* set macrophage type */
  set_type(tp);

  MacrophageState state;
  if (st == MacrophageState::null) {
    /* calculate random number between 0 and 1 */
    double r = distribution_0_1(rng);
    if (r < params.macrophage_active_probability) {
      /* macrophage state is active */
      state = (validation_chemotaxis_and_recruitment) ? MacrophageState::active : MacrophageState::resting;
    } else {
      /* macrophage state is resting */
      state = MacrophageState::resting;
    }
  } else {
    state = st;
  }

  total_virus_uptaken = 0.0;

  /* resident on cell */
  resident_on_cell = nullptr;

  /* set macrophage state, age, lifespan, secretions and uptake */
  set_state(state);

  /* set the area */
  double mac_area = CAM_PI * params.macrophage_radius * params.macrophage_radius;
  set_area(mac_area);

  /* update the reaction diffusion secretions & uptakes */
  inc_rds_values();

  /* apoptosis */
  is_being_apoptosed = false;
  apoptosis_start_time = 0.0;
  apoptosis_timescale = params.macrophage_apoptosis_timescale;
  num_cells_apoptosing = 0;
  is_apoptosing_epithelial = false;
  apoptosing_agent = nullptr;

  /* phagocytosis */
  phagocytosis_timescale = params.macrophage_phagocytosis_timescale;
  is_phagocytosing_epithelial = false;
  phagocytosing_agent = nullptr;
  phagocytosed_by = nullptr;

  std::string outf;
  if (use_output) {
    outf = output_folder + std::string("/path_macrophage.txt");
    output = fopen(outf.c_str(), "w");

    outf = output_folder + std::string("/macrophage_total_virus_uptaken.txt");
    uptake_output = fopen(outf.c_str(), "w");

    outf = output_folder + std::string("/macrophage_lifespan.txt");
    lifespan_output = fopen(outf.c_str(), "w");

    outf = output_folder + std::string("/macrophage_movement_rate.txt");
    mrate_output = fopen(outf.c_str(), "w");
  } else {
    output = nullptr;
    uptake_output = nullptr;
    lifespan_output = nullptr;
    mrate_output = nullptr;
  }

  /* successful */
  (*status) = CAM_SUCCESS;
}


/*
 * =============================================================================
 * Macrophage destructor
 * =============================================================================
 */
Macrophage::~Macrophage ()
{
  /* need to decrease the cytokine secretions and uptakes */
  dec_rds_values();

  if (resident_on_cell != nullptr) {
    resident_on_cell = nullptr;
  }

  if (output != nullptr) {
    fclose(output);
  }

  if (uptake_output != nullptr) {
    fclose(uptake_output);
  }

  if (lifespan_output != nullptr) {
    fclose(lifespan_output);
  }

  if (mrate_output != nullptr) {
    fclose(mrate_output);
  }
}


/*
 * =============================================================================
 * set macrophage state, reset age, lifespan, secretions and uptakes
 * =============================================================================
 */
void
Macrophage::set_state (const MacrophageState st)
{
  /* set state */
  mac_state = st;

  /* reset age */
  set_age(0.0);

  /* set the lifespan, depending on type and state */
  double lspan = 0.0;
  switch (mac_state) {
    case MacrophageState::resting:
      lspan = params.macrophage_resting_lifespan;
      break;
    case MacrophageState::active:
      lspan = params.macrophage_active_lifespan;
      break;
    default:
      break;
  }
  set_lifespan(lspan, true);

  /* set secretion depending on type and state */
  set_secretion(st);

  /* set uptake depending on type and state */
  set_uptake(st);

  /* set the movement rate & bias, depending on state and type */
  set_movement(st);
}


/*
 * =============================================================================
 * set the cytokine and chemokine secretion based on state
 * =============================================================================
 */
void
Macrophage::set_secretion (const MacrophageState st)
{
  switch (st) {
    case MacrophageState::resting:
      basal_secreted_ifn_1 = params.macrophage_resting_basal_secreted_ifn_1;
      secreted_ifn_1 = params.macrophage_resting_secreted_ifn_1;

      basal_secreted_il_6 = params.macrophage_resting_basal_secreted_il_6;
      secreted_il_6 = params.macrophage_resting_secreted_il_6;

      basal_secreted_il_10 = params.macrophage_resting_basal_secreted_il_10;
      secreted_il_10 = params.macrophage_resting_secreted_il_10;

      basal_secreted_il_12 = params.macrophage_resting_basal_secreted_il_12;
      secreted_il_12 = params.macrophage_resting_secreted_il_12;

      basal_secreted_g_csf = params.macrophage_resting_basal_secreted_g_csf;
      secreted_g_csf = params.macrophage_resting_secreted_g_csf;

      basal_secreted_m_csf = params.macrophage_resting_basal_secreted_m_csf;
      secreted_m_csf = params.macrophage_resting_secreted_m_csf;

      basal_secreted_epi_inf_ck = 0.0;
      secreted_epi_inf_ck = 0.0;

      secreted_apop_ck = params.macrophage_secreted_apop_ck;
      secreted_phago_ck = 0.0;

      break;
    case MacrophageState::active:
      basal_secreted_ifn_1 = params.macrophage_active_basal_secreted_ifn_1;
      secreted_ifn_1 = params.macrophage_active_secreted_ifn_1;

      basal_secreted_il_6 = params.macrophage_active_basal_secreted_il_6;
      secreted_il_6 = params.macrophage_active_secreted_il_6;

      basal_secreted_il_10 = params.macrophage_active_basal_secreted_il_10;
      secreted_il_10 = params.macrophage_active_secreted_il_10;

      basal_secreted_il_12 = params.macrophage_active_basal_secreted_il_12;
      secreted_il_12 = params.macrophage_active_secreted_il_12;

      basal_secreted_g_csf = params.macrophage_active_basal_secreted_g_csf;
      secreted_g_csf = params.macrophage_active_secreted_g_csf;

      basal_secreted_m_csf = params.macrophage_active_basal_secreted_m_csf;
      secreted_m_csf = params.macrophage_active_secreted_m_csf;

      basal_secreted_epi_inf_ck = 0.0;
      secreted_epi_inf_ck = 0.0;

      secreted_apop_ck = params.macrophage_secreted_apop_ck;
      secreted_phago_ck = 0.0;

      break;
    case MacrophageState::apoptotic:
      basal_secreted_ifn_1 = 0.0;
      secreted_ifn_1 = 0.0;

      basal_secreted_il_6 = 0.0;
      secreted_il_6 = 0.0;

      basal_secreted_il_10 = 0.0;
      secreted_il_10 = 0.0;

      basal_secreted_il_12 = 0.0;
      secreted_il_12 = 0.0;

      basal_secreted_g_csf = 0.0;
      secreted_g_csf = 0.0;

      basal_secreted_m_csf = 0.0;
      secreted_m_csf = 0.0;

      basal_secreted_epi_inf_ck = 0.0;
      secreted_epi_inf_ck = 0.0;

      secreted_apop_ck = 0.0;
      secreted_phago_ck = params.macrophage_apoptotic_secreted_phago_ck;

      break;
    default:
      basal_secreted_ifn_1 = 0.0;
      secreted_ifn_1 = 0.0;

      basal_secreted_il_6 = 0.0;
      secreted_il_6 = 0.0;

      basal_secreted_il_10 = 0.0;
      secreted_il_10 = 0.0;

      basal_secreted_il_12 = 0.0;
      secreted_il_12 = 0.0;

      basal_secreted_g_csf = 0.0;
      secreted_g_csf = 0.0;

      basal_secreted_m_csf = 0.0;
      secreted_m_csf = 0.0;

      basal_secreted_epi_inf_ck = 0.0;
      secreted_epi_inf_ck = 0.0;

      secreted_apop_ck = 0.0;
      secreted_phago_ck = 0.0;

      break;
  }

  regulate_secretion();
}


/*
 * =============================================================================
 * set the cytokine and chemokine uptake based on state
 * =============================================================================
 */
void
Macrophage::set_uptake (const MacrophageState st)
{
  switch (st) {
    case MacrophageState::resting:
      uptake_ifn_1 = params.macrophage_resting_uptake_ifn_1;

      uptake_il_6 = params.macrophage_resting_uptake_il_6;
      uptake_il_10 = params.macrophage_resting_uptake_il_10;
      uptake_il_12 = params.macrophage_resting_uptake_il_12;

      uptake_g_csf = params.macrophage_resting_uptake_g_csf;
      uptake_m_csf = params.macrophage_resting_uptake_m_csf;

      uptake_epi_inf_ck = params.macrophage_resting_uptake_epi_inf_ck;
      uptake_apop_ck = 0.0;
      uptake_phago_ck = params.macrophage_resting_uptake_phago_ck;

      uptake_virus_inf = params.macrophage_resting_uptake_virus_inf;
      decay_virus_inf = params.macrophage_resting_decay_virus_inf;

      break;
    case MacrophageState::active:
      uptake_ifn_1 = params.macrophage_active_uptake_ifn_1;

      uptake_il_6 = params.macrophage_active_uptake_il_6;
      uptake_il_10 = params.macrophage_active_uptake_il_10;
      uptake_il_12 = params.macrophage_active_uptake_il_12;

      uptake_g_csf = params.macrophage_active_uptake_g_csf;
      uptake_m_csf = params.macrophage_active_uptake_m_csf;

      uptake_epi_inf_ck = params.macrophage_active_uptake_epi_inf_ck;
      uptake_apop_ck = 0.0;
      uptake_phago_ck = params.macrophage_active_uptake_phago_ck;

      uptake_virus_inf = params.macrophage_active_uptake_virus_inf;
      decay_virus_inf = params.macrophage_active_decay_virus_inf;

      break;
    default:
      uptake_ifn_1 = 0.0;

      uptake_il_6 = 0.0;
      uptake_il_10 = 0.0;
      uptake_il_12 = 0.0;

      uptake_g_csf = 0.0;
      uptake_m_csf = 0.0;

      uptake_epi_inf_ck = 0.0;
      uptake_apop_ck = 0.0;
      uptake_phago_ck = 0.0;

      uptake_virus_inf = 0.0;
      decay_virus_inf = 0.0;

      break;
  }

  regulate_uptake();
}


void
Macrophage::set_movement (const MacrophageState st)
{
  switch (st) {
    case MacrophageState::resting:
      movement_rate = params.macrophage_resting_movement_rate;
      movement_bias = params.macrophage_resting_movement_bias;
      break;
    case MacrophageState::active:
      movement_rate = params.macrophage_active_movement_rate;
      movement_bias = params.macrophage_active_movement_bias;
      break;
    default:
      movement_rate = 0.0;
      movement_bias = 0.0;
      break;
  }
}


void
Macrophage::change_state (const MacrophageState st)
{
  /* decrease the secretions and uptakes */
  dec_rds_values();

  /* set the new state */
  /* also resets the age, cytokines secretions & uptakes, movement rate and bias */
  set_state(st);

  /* increase the secretions and uptakes */
  inc_rds_values();
}


void
Macrophage::update_lifespan (const MacrophageState st)
{
  (void) st;

  return; // currently not used

  double inhib = cam_util_hill_function(
    params.macrophage_lifespan_downreg_virus_EC50, total_virus_uptaken, params.macrophage_lifespan_downreg_virus_hill
  );

  double lspan = get_lifespan(true)*inhib; // get_lifespan(true) returns max lifespan
  set_lifespan(lspan, false);
}


void
Macrophage::update_movement (const MacrophageState st)
{
  /* reset the movement rate */
  set_movement(st);

  return; // currently not used

  double conc = cam_util_hill_function(
    total_virus_uptaken, params.macrophage_movement_upreg_virus_EC50, params.macrophage_movement_upreg_virus_hill
  );

  movement_rate *= (1.0 + conc); // NEEDSWORK: => max at 2*movement_rate
}


void
Macrophage::regulate_secretion ()
{
  double conc = 0.0;
  double signal = 0.0;
  double upreg = 0.0;
  double downreg = 0.0;

  std::pair<int, int> xy = get_coordinates();

  const double ifn_1 = rd_sys->ifn_1->get_conc(xy.first, xy.second);
  const double il_6 = rd_sys->il_6->get_conc(xy.first, xy.second);
  const double il_10 = rd_sys->il_10->get_conc(xy.first, xy.second);
  const double il_12 = rd_sys->il_12->get_conc(xy.first, xy.second);
  const double g_csf = rd_sys->g_csf->get_conc(xy.first, xy.second);
  const double m_csf = rd_sys->m_csf->get_conc(xy.first, xy.second);

  /* T1 IFN */
  signal = cam_util_hill_function(
    total_virus_uptaken, params.macrophage_ifn_1_signal_virus_EC50, params.macrophage_ifn_1_signal_virus_hill
  );

  downreg = cam_util_hill_function(
    params.macrophage_ifn_1_downreg_il_10_EC50, il_10, params.macrophage_ifn_1_downreg_il_10_hill
  );

  conc = ifn_1;
  upreg = cam_util_hill_function(
    ifn_1, params.macrophage_ifn_1_EC50, params.macrophage_ifn_1_hill
  );

  secreted_ifn_1 = (
    ((basal_secreted_ifn_1 / (1.0 + conc*conc)) + secreted_ifn_1*upreg*downreg)*signal
  );


  /* IL-6 */
  signal = cam_util_hill_function(
    total_virus_uptaken, params.macrophage_il_6_signal_virus_EC50, params.macrophage_il_6_signal_virus_hill
  );

  downreg = cam_util_hill_function(
    params.macrophage_il_6_downreg_il_10_EC50, il_10, params.macrophage_il_6_downreg_il_10_hill
  );

  conc = il_6;
  upreg = cam_util_hill_function(
    il_6, params.macrophage_il_6_EC50, params.macrophage_il_6_hill
  );

  secreted_il_6 = (
    ((basal_secreted_il_6 / (1.0 + conc*conc)) + secreted_il_6*upreg*downreg)*signal
  );


  /* IL-12 */
  signal = cam_util_hill_function(
    total_virus_uptaken, params.macrophage_il_12_signal_virus_EC50, params.macrophage_il_12_signal_virus_hill
  );

  downreg = cam_util_hill_function(
    params.macrophage_il_12_downreg_il_10_EC50, il_10, params.macrophage_il_12_downreg_il_10_hill
  );

  conc = il_12;
  upreg = cam_util_hill_function(
    il_12, params.macrophage_il_12_EC50, params.macrophage_il_12_hill
  );

  secreted_il_12 = (
    ((basal_secreted_il_12 / (1.0 + conc*conc)) + secreted_il_12*upreg*downreg)*signal
  );


  /* IL-10 */
  signal = 0.25 * cam_util_hill_function(
    il_6, params.macrophage_il_10_signal_il_6_EC50, params.macrophage_il_10_signal_il_6_hill
  );
  signal += 0.25 * cam_util_hill_function(
    il_12, params.macrophage_il_10_signal_il_12_EC50, params.macrophage_il_10_signal_il_12_hill
  );
  signal += 0.25 * cam_util_hill_function(
    g_csf, params.macrophage_il_10_signal_g_csf_EC50, params.macrophage_il_10_signal_g_csf_hill
  );
  signal += 0.25 * cam_util_hill_function(
    m_csf, params.macrophage_il_10_signal_m_csf_EC50, params.macrophage_il_10_signal_m_csf_hill
  );

  conc = il_10;
  upreg = cam_util_hill_function(
    il_10, params.macrophage_il_10_EC50, params.macrophage_il_10_hill
  );

  secreted_il_10 = (
    ((basal_secreted_il_10 / (1.0 + conc*conc)) + secreted_il_10*upreg)*signal
  );


  /* G-CSF */
  signal = cam_util_hill_function(
    total_virus_uptaken, params.macrophage_g_csf_signal_virus_EC50, params.macrophage_g_csf_signal_virus_hill
  );

  downreg = cam_util_hill_function(
    params.macrophage_g_csf_downreg_il_10_EC50, il_10, params.macrophage_g_csf_downreg_il_10_hill
  );

  conc = g_csf;
  upreg = cam_util_hill_function(
    g_csf, params.macrophage_g_csf_EC50, params.macrophage_g_csf_hill
  );

  secreted_g_csf = (
    ((basal_secreted_g_csf / (1.0 + conc*conc)) + secreted_g_csf*upreg*downreg)*signal
  );


  /* M-CSF */
  signal = cam_util_hill_function(
    total_virus_uptaken, params.macrophage_m_csf_signal_virus_EC50, params.macrophage_m_csf_signal_virus_hill
  );

  downreg = cam_util_hill_function(
    params.macrophage_m_csf_downreg_il_10_EC50, il_10, params.macrophage_m_csf_downreg_il_10_hill
  );

  conc = m_csf;
  upreg = cam_util_hill_function(
    m_csf, params.macrophage_m_csf_EC50, params.macrophage_m_csf_hill
  );

  secreted_m_csf = (
    ((basal_secreted_m_csf / (1.0 + conc*conc)) + secreted_m_csf*upreg*downreg)*signal
  );


  /* infected epithelial chemokine - not used */
}


void
Macrophage::regulate_uptake ()
{
  std::pair<int, int> xy = get_coordinates();

  double epi_area = ( (resident_on_cell == nullptr) ? 0.0 : resident_on_cell->get_max_area());
  double n_v = rd_sys->virus_inf->get_conc(xy.first, xy.second) * epi_area;

  double gfnc = 0.0;
  if (params.macrophage_blocking_function) {
    gfnc = 1.0;
  } else if (n_v >= 1.0) {
    params.macrophage_blocking_function = true;
    gfnc = 1.0;
  }

  uptake_virus_inf *= n_v * gfnc * cam_util_hill_function(
    params.macrophage_virus_inf_downreg_internalised_virus_EC50,
    total_virus_uptaken,
    params.macrophage_virus_inf_downreg_internalised_virus_hill
  );
}


void
Macrophage::update_secretions ()
{
  /* update the reaction diffusion secretions and uptakes */
  dec_rds_values();

  /* update the lifespan */
  update_lifespan(mac_state);

  /* reset the secretions */
  set_secretion(mac_state);

  /* reset the uptakes */
  set_uptake(mac_state);

  /* update the movement */
  update_movement(mac_state);

  /* update the reaction diffusion secretions and uptakes */
  inc_rds_values();
}


/*
 * =============================================================================
 * Compute the macrophage action (activation, deactivation, phagocytosis)
 * =============================================================================
 */
void
Macrophage::action (const double time)
{
  /* random number */
  double r;

  /* coords */
  std::pair<int, int> xy = get_coordinates();

  if (mac_state == MacrophageState::apoptotic || is_being_apoptosed || is_being_phagocytosed || mac_state == MacrophageState::phagocytosed) {
    return;
  }

  /* update the secretions */
  update_secretions();

  double ifn_uptake = rd_sys->ifn_1->get_conc(xy.first, xy.second);
  double trigger_apop = params.macrophage_apoptosis_by_ifn_max_rate * cam_util_hill_function(ifn_uptake, params.macrophage_apoptosis_by_ifn_half_max, params.macrophage_apoptosis_by_ifn_hill);

  r = distribution_0_1(rng);
  if (r < trigger_apop) {
    increase_age(10.0*get_lifespan(true));

    return;
  }

  /* advance internalised virus content */
  total_virus_uptaken += timestep * (uptake_virus_inf - decay_virus_inf * total_virus_uptaken); // virions
  if (params.macrophage_activate_weight_virus_inf * total_virus_uptaken >= 1.0 && mac_state == MacrophageState::resting) {
    /* if total_virus_uptaken is beyond the threshold then burst */
    if (total_virus_uptaken > params.macrophage_resting_burst_threshold) {
      /* flag for bursting - handled in environment */
      flag_bursting = true;
      return;
    }

    /* activation */
    r = distribution_0_1(rng);
    if (r < params.macrophage_active_probability && !is_phagocytosing_epithelial && phagocytosing_agent == nullptr && !is_being_apoptosed && !is_being_phagocytosed) {
      change_state(MacrophageState::active);

      /* update the secretions */
      update_secretions();

      return;
    }
  }


  /* get the aggregate of the cytokines */
  double cyt_tmp = (
    params.macrophage_activate_weight_ifn_1 * rd_sys->ifn_1->get_conc(xy.first, xy.second) +
    params.macrophage_activate_weight_il_6 * rd_sys->il_6->get_conc(xy.first, xy.second) +
    params.macrophage_activate_weight_m_csf * rd_sys->m_csf->get_conc(xy.first, xy.second) +
    params.macrophage_activate_weight_epi_inf_ck * rd_sys->epi_inf_ck->get_conc(xy.first, xy.second) +
    params.macrophage_activate_weight_phago_ck * rd_sys->phago_ck->get_conc(xy.first, xy.second) -
    params.macrophage_activate_weight_il_10 * rd_sys->il_10->get_conc(xy.first, xy.second)
  );

  double cyt = cam_util_hill_function(
    cyt_tmp, params.macrophage_activate_EC50, params.macrophage_activate_hill
  );

  /* resting actions */
  if (mac_state == MacrophageState::resting) {
    /* if total_virus_uptaken is beyond the threshold then burst */
    if (total_virus_uptaken > params.macrophage_resting_burst_threshold) {
      /* flag for bursting - handled in environment */
      flag_bursting = true;
      return;
    }

    /* activation */
    r = distribution_0_1(rng);
    if (r < cyt && !is_phagocytosing_epithelial && phagocytosing_agent == nullptr && !is_being_apoptosed && !is_being_phagocytosed) {
      change_state(MacrophageState::active);

      /* update the secretions */
      update_secretions();

      return;
    }
  }


  /* active actions */
  if (mac_state == MacrophageState::active) {
    /* if total_virus_uptaken is beyond the threshold then burst */
    if (total_virus_uptaken > params.macrophage_active_burst_threshold) {
      /* flag for bursting - handled in environment */
      flag_bursting = true;
      return;
    }

    /* phagocytosis of epithelial */
    int ret_code = phagocytose_epithelial(time);
    if (ret_code == 0) {
      /* successful phagocytosis action; start, or end */
      return;
    }

    /* phagocytosis of apoptotic immune cells */
    phagocytose_immune(time);

    /* deactivation */
    r = distribution_0_1(rng);
    double r2 = distribution_0_1(rng);
    if (r > cyt && r2 > params.macrophage_active_probability && params.macrophage_activate_weight_virus_inf * total_virus_uptaken < 1.0 &&
        !is_phagocytosing_epithelial && phagocytosing_agent == nullptr && !is_being_apoptosed && !is_being_phagocytosed) {
      change_state(MacrophageState::resting);

      /* update the secretions */
      update_secretions();

      return;
    }
  }
}


int
Macrophage::phagocytose_epithelial (const double time)
{
  if (phagocytosing_agent != nullptr) {
    /* this macrophage cell is already phagocytosing another immune cell */
    return -1;
  }


  /* phagocytosis of epithelial cell */
  if (is_phagocytosing_epithelial && !resident_on_cell->is_being_phagocytosed) {
    /* phagocytosis is complete */

    is_phagocytosing_epithelial = false;

    return 0;

  } else if (!is_phagocytosing_epithelial && !resident_on_cell->is_being_phagocytosed) {
    /* macrophage starting phagocytosis */

    if (resident_on_cell->check_state(EpithelialState::apoptotic)) {
      double r = distribution_0_1(rng);
      if (r < params.macrophage_active_phagocytose_apoptotic_epithelial_probability) {
        /* start phagocytosis */

        double tscale = phagocytosis_timescale;

        resident_on_cell->start_phagocytosis(time);
        resident_on_cell->set_phagocytosis_timescale(tscale);
        resident_on_cell->switch_flag_phagocytosis();
        is_phagocytosing_epithelial = true;

        return 0;
      }
    } else if (resident_on_cell->check_state(EpithelialState::infected) || resident_on_cell->check_state(EpithelialState::eclipse)) {
      double r = distribution_0_1(rng);
      if (r < params.macrophage_active_phagocytose_infected_epithelial_probability) {
        /* start phagocytosis */

        double tscale = phagocytosis_timescale;

        resident_on_cell->start_phagocytosis(time);
        resident_on_cell->set_phagocytosis_timescale(tscale);
        resident_on_cell->switch_flag_phagocytosis();
        is_phagocytosing_epithelial = true;

        return 0;
      }
    }
  }

  /* want to check whether immune cell apoptosis can be triggered */
  return -2;
}


void
Macrophage::phagocytose_immune (const double time)
{
  if (is_phagocytosing_epithelial) {
    /* this macrophage is phagocytosing an epithelial cell then return */
    return;
  }


  /* phagocytosis of resident immune cell */
  if (phagocytosing_agent == nullptr) {
    /* starting phagocytosis of immune cell */

    std::vector<Agent*>* agents = resident_on_cell->get_residents();
    for (int i = 0; i < (int) (*agents).size(); i++) {
      Agent* ag = (*agents)[i];
      if (ag == this) {
        /* NOTE: might have to dynamic_cast "this" object to Agent* */
        /* if the agent is the current object */
        continue;
      }

      if (!ag->is_apoptotic()) {
        /* only apoptotic immune cells can be phagocytosed */
        continue;
      }

      double r = distribution_0_1(rng);
      if (r < (params.macrophage_active_phagocytose_immune_probability)) {
        /* start phagocytosis */

        double tscale = phagocytosis_timescale;

        ag->start_phagocytosis(time);
        ag->set_phagocytosis_timescale(tscale);

        phagocytosing_agent = ag;
        ag->phagocytosed_by = this;

        return;
      }
    }

  } else {
    if (!phagocytosing_agent->is_being_phagocytosed) {
      /* end of phagocytosis */

      phagocytosing_agent = nullptr;
    } /* else { phagocytosis ongoing } */
  }
}


int
Macrophage::migrate (const double dt, const double tm, std::function<int(Agent*, std::vector<std::pair<int,int>>)> fndsp, std::function<bool(Agent*, std::pair<int,int>)> chksp)
{
  Epithelial* epi = resident_on_cell;
  const int sz = (int) epi->immediate_moore.size();

  /* return index */
  int ind = sz;

  /* movement rate */
  if (is_phagocytosing_epithelial || phagocytosing_agent != nullptr || is_being_phagocytosed ||
      is_being_apoptosed || (mac_state == MacrophageState::apoptotic)) {
    return ind;
  }

  if (movement_rate < CAM_TOL) {
    return ind;
  }

  /* movement */
  if (((int) round(tm/dt) % (int) round(movement_rate/dt)) == 0) {
    std::pair<int, int> xy;

    /* Get the cytokine levels from all immediate neighbours */
    double total_cyt = 0.0;
    double nbr_cyt = 0.0;

    for (int j = 0; j < sz; j++) {
      xy = epi->immediate_moore[j];

      /* cytokines/chemokines used for migration */
      nbr_cyt = (
        params.macrophage_movement_weight_epi_inf_ck * rd_sys->epi_inf_ck->get_conc(xy.first, xy.second) +
        params.macrophage_movement_weight_phago_ck * rd_sys->phago_ck->get_conc(xy.first, xy.second)
      );

      total_cyt += pow(nbr_cyt, movement_bias);
    }

    if (total_cyt < CAM_TOL || mac_state == MacrophageState::resting) {
      /* RANDOM MIGRATION: choose an immediate neighbour */

      /* find space using function provided as arg */
      ind = fndsp(this, epi->immediate_moore);
      if (ind < 0)
        ind = sz;
    } else {
      /* CHEMOTACTIC MIGRATION: (biased random walk) choose neighbour based on cytokine levels */

      double r = distribution_0_1(rng)*total_cyt;
      double running_total = 0.0;
      for (int j = 0; j < sz; j++) {
        /* coordinates of neighbour */
        xy = epi->immediate_moore[j];

        /* cytokine/chemokine concentration used for migration */
        nbr_cyt = (
          params.macrophage_movement_weight_epi_inf_ck * rd_sys->epi_inf_ck->get_conc(xy.first, xy.second) +
          params.macrophage_movement_weight_phago_ck * rd_sys->phago_ck->get_conc(xy.first, xy.second)
        );

        running_total += pow(nbr_cyt, movement_bias);

        bool successful = chksp(this, xy);
        if (running_total >= r && successful) {
          ind = j;
          break;
        }
      }
    }
  }

  return ind;
}


bool
Macrophage::recruit (std::vector<std::pair<int,int>> nbrs)
{
  if (params.macrophage_recruitment_hill < CAM_TOL) {
    return false;
  }

  /* Get the cytokine levels from all neighbours */
  double total_cyt = 0.0;
  double nbr_m_csf = 0.0, nbr_il_6 = 0.0, nbr_epi_inf_ck = 0.0, nbr_phago_ck = 0.0, nbr_il_10 = 0.0;

  for (std::pair<int,int> nb:nbrs) {
    nbr_m_csf += rd_sys->m_csf->get_conc(nb.first, nb.second);
    nbr_il_6 += rd_sys->il_6->get_conc(nb.first, nb.second);
    nbr_epi_inf_ck += rd_sys->epi_inf_ck->get_conc(nb.first, nb.second);
    nbr_phago_ck += rd_sys->phago_ck->get_conc(nb.first, nb.second);
    nbr_il_10 += rd_sys->il_10->get_conc(nb.first, nb.second);
  }

  total_cyt = (
    params.macrophage_recruitment_weight_m_csf * nbr_m_csf +
    params.macrophage_recruitment_weight_il_6 * nbr_il_6 +
    params.macrophage_recruitment_weight_epi_inf_ck * nbr_epi_inf_ck +
    params.macrophage_recruitment_weight_phago_ck * nbr_phago_ck -
    params.macrophage_recruitment_weight_il_10 * nbr_il_10
  ) / ((double) nbrs.size());

  double cyt = cam_util_hill_function(
    total_cyt, params.macrophage_recruitment_EC50, params.macrophage_recruitment_hill
  );

  double r = distribution_0_1(rng);
  if (r < cyt) {
    return true;
  }

  return false;
}


/*
 * =============================================================================
 * Start the apoptosis in the macrophage
 * =============================================================================
 */
void
Macrophage::start_apoptosis (const double time)
{
  std::pair<int, int> xy = get_coordinates();

  is_being_apoptosed = true;
  apoptosis_start_time = time;
  num_cells_apoptosing = 1;

  /* secrete the chemokine which recruits NK cells to aid in apoptosis */
  rd_sys->apop_ck->inc_source(xy.first, xy.second, secreted_apop_ck, true);
}


/*
 * =============================================================================
 * End the apoptosis in the macrophage
 * =============================================================================
 */
void
Macrophage::end_apoptosis ()
{
  std::pair<int, int> xy = get_coordinates();

  is_being_apoptosed = false;
  apoptosis_start_time = 0.0;
  num_cells_apoptosing = 0;

  change_state(MacrophageState::apoptotic);

  /* update the secretions */
  update_secretions();

  /* decrease the NK chemokine */
  rd_sys->apop_ck->dec_source(xy.first, xy.second, secreted_apop_ck, true);
}


/*
 * =============================================================================
 * Set the apoptosis timescale
 * =============================================================================
 */
void
Macrophage::set_apoptosis_timescale (const double tms)
{
  double timescale = tms / ((double) num_cells_apoptosing);

  /* change time scale of apoptosis due to NK action */
  if (timescale < apoptosis_timescale) { 
    apoptosis_timescale = timescale;
  }
}


/*
 * =============================================================================
 * Start the phagocytosis in the macrophage
 * =============================================================================
 */
void
Macrophage::start_phagocytosis (const double time)
{
  is_being_phagocytosed = true;
  phagocytosis_start_time = time;
}


/*
 * =============================================================================
 * End the phagocytosis in the macrophage
 * =============================================================================
 */
void
Macrophage::end_phagocytosis ()
{
  is_being_phagocytosed = false;
  phagocytosis_start_time = 0.0;

  phagocytosed_by->phagocytosing_agent = nullptr;
  phagocytosed_by = nullptr;
  phagocytosing_agent = nullptr;

  change_state(MacrophageState::phagocytosed);

  /* update the secretions */
  update_secretions();
}


/*
 * =============================================================================
 * Set the phagocytosis timescale
 * =============================================================================
 */
void
Macrophage::set_phagocytosis_timescale (const double tms)
{
  (void) tms; // phagocytosis is set by the target cell only

  return;

  /* change time scale of phagocytosis */
  if (tms < phagocytosis_timescale) {
    phagocytosis_timescale = tms;
  }
}
