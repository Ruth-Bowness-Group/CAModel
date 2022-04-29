#include <cstring>

#include "cam_monocyte.h"
#include "monocyte.h"
#include "../../cam_error.h"
#include "../../cam_random.h"
#include "../../cam_util.h"
#include "../../cam_parameters.h"


/* struct containing monocyte specific parameters */
static monocyte_params params;


/*
 * =============================================================================
 * Set the monocyte specific parameters
 * =============================================================================
 */
int
cam_monocyte_set_params (YAML::Node* cfg)
{
  // int status = 0;
  std::string ren;

  try {
    /* radius */
    params.monocyte_radius = (*cfg)["monocyte_radius"].as<double>();

    /* probability of activating monocyte */
    params.monocyte_active_probability = (*cfg)["monocyte_active_probability"].as<double>();

    /* lifespan */
    params.monocyte_resting_lifespan = (*cfg)["monocyte_resting_lifespan"].as<double>();
    params.monocyte_active_lifespan = (*cfg)["monocyte_active_lifespan"].as<double>();

    params.monocyte_lifespan_downreg_virus_EC50 = (*cfg)["monocyte_lifespan_downreg_virus_EC50"].as<double>();
    params.monocyte_lifespan_downreg_virus_hill = (*cfg)["monocyte_lifespan_downreg_virus_hill"].as<double>();

    /*
     * cytokines
     */

    /* T1 IFN */
    params.monocyte_resting_basal_secreted_ifn_1 = (*cfg)["monocyte_resting_basal_secreted_ifn_1"].as<double>();
    params.monocyte_resting_secreted_ifn_1 = (*cfg)["monocyte_resting_secreted_ifn_1"].as<double>();
    params.monocyte_resting_uptake_ifn_1 = (*cfg)["monocyte_resting_uptake_ifn_1"].as<double>();

    params.monocyte_active_basal_secreted_ifn_1 = (*cfg)["monocyte_active_basal_secreted_ifn_1"].as<double>();
    params.monocyte_active_secreted_ifn_1 = (*cfg)["monocyte_active_secreted_ifn_1"].as<double>();
    params.monocyte_active_uptake_ifn_1 = (*cfg)["monocyte_active_uptake_ifn_1"].as<double>();

    params.monocyte_ifn_1_EC50 = (*cfg)["monocyte_ifn_1_EC50"].as<double>();
    params.monocyte_ifn_1_hill = (*cfg)["monocyte_ifn_1_hill"].as<double>();
    params.monocyte_ifn_1_signal_virus_EC50 = (*cfg)["monocyte_ifn_1_signal_virus_EC50"].as<double>();
    params.monocyte_ifn_1_signal_virus_hill = (*cfg)["monocyte_ifn_1_signal_virus_hill"].as<double>();
    params.monocyte_ifn_1_downreg_il_10_EC50 = (*cfg)["monocyte_ifn_1_downreg_il_10_EC50"].as<double>();
    params.monocyte_ifn_1_downreg_il_10_hill = (*cfg)["monocyte_ifn_1_downreg_il_10_hill"].as<double>();

    /* IL-6 */
    params.monocyte_resting_basal_secreted_il_6 = (*cfg)["monocyte_resting_basal_secreted_il_6"].as<double>();
    params.monocyte_resting_secreted_il_6 = (*cfg)["monocyte_resting_secreted_il_6"].as<double>();
    params.monocyte_resting_uptake_il_6 = (*cfg)["monocyte_resting_uptake_il_6"].as<double>();

    params.monocyte_active_basal_secreted_il_6 = (*cfg)["monocyte_active_basal_secreted_il_6"].as<double>();
    params.monocyte_active_secreted_il_6 = (*cfg)["monocyte_active_secreted_il_6"].as<double>();
    params.monocyte_active_uptake_il_6 = (*cfg)["monocyte_active_uptake_il_6"].as<double>();

    params.monocyte_il_6_EC50 = (*cfg)["monocyte_il_6_EC50"].as<double>();
    params.monocyte_il_6_hill = (*cfg)["monocyte_il_6_hill"].as<double>();
    params.monocyte_il_6_signal_virus_EC50 = (*cfg)["monocyte_il_6_signal_virus_EC50"].as<double>();
    params.monocyte_il_6_signal_virus_hill = (*cfg)["monocyte_il_6_signal_virus_hill"].as<double>();
    params.monocyte_il_6_downreg_il_10_EC50 = (*cfg)["monocyte_il_6_downreg_il_10_EC50"].as<double>();
    params.monocyte_il_6_downreg_il_10_hill = (*cfg)["monocyte_il_6_downreg_il_10_hill"].as<double>();

    /* IL-10 */
    params.monocyte_resting_basal_secreted_il_10 = (*cfg)["monocyte_resting_basal_secreted_il_10"].as<double>();
    params.monocyte_resting_secreted_il_10 = (*cfg)["monocyte_resting_secreted_il_10"].as<double>();
    params.monocyte_resting_uptake_il_10 = (*cfg)["monocyte_resting_uptake_il_10"].as<double>();

    params.monocyte_active_basal_secreted_il_10 = (*cfg)["monocyte_active_basal_secreted_il_10"].as<double>();
    params.monocyte_active_secreted_il_10 = (*cfg)["monocyte_active_secreted_il_10"].as<double>();
    params.monocyte_active_uptake_il_10 = (*cfg)["monocyte_active_uptake_il_10"].as<double>();

    params.monocyte_il_10_EC50 = (*cfg)["monocyte_il_10_EC50"].as<double>();
    params.monocyte_il_10_hill = (*cfg)["monocyte_il_10_hill"].as<double>();
    params.monocyte_il_10_signal_il_6_EC50 = (*cfg)["monocyte_il_10_signal_il_6_EC50"].as<double>();
    params.monocyte_il_10_signal_il_6_hill = (*cfg)["monocyte_il_10_signal_il_6_hill"].as<double>();
    params.monocyte_il_10_signal_il_12_EC50 = (*cfg)["monocyte_il_10_signal_il_12_EC50"].as<double>();
    params.monocyte_il_10_signal_il_12_hill = (*cfg)["monocyte_il_10_signal_il_12_hill"].as<double>();
    params.monocyte_il_10_signal_g_csf_EC50 = (*cfg)["monocyte_il_10_signal_g_csf_EC50"].as<double>();
    params.monocyte_il_10_signal_g_csf_hill = (*cfg)["monocyte_il_10_signal_g_csf_hill"].as<double>();
    params.monocyte_il_10_signal_m_csf_EC50 = (*cfg)["monocyte_il_10_signal_m_csf_EC50"].as<double>();
    params.monocyte_il_10_signal_m_csf_hill = (*cfg)["monocyte_il_10_signal_m_csf_hill"].as<double>();

    /* IL-12 */
    params.monocyte_resting_basal_secreted_il_12 = (*cfg)["monocyte_resting_basal_secreted_il_12"].as<double>();
    params.monocyte_resting_secreted_il_12 = (*cfg)["monocyte_resting_secreted_il_12"].as<double>();
    params.monocyte_resting_uptake_il_12 = (*cfg)["monocyte_resting_uptake_il_12"].as<double>();

    params.monocyte_active_basal_secreted_il_12 = (*cfg)["monocyte_active_basal_secreted_il_12"].as<double>();
    params.monocyte_active_secreted_il_12 = (*cfg)["monocyte_active_secreted_il_12"].as<double>();
    params.monocyte_active_uptake_il_12 = (*cfg)["monocyte_active_uptake_il_12"].as<double>();

    params.monocyte_il_12_EC50 = (*cfg)["monocyte_il_12_EC50"].as<double>();
    params.monocyte_il_12_hill = (*cfg)["monocyte_il_12_hill"].as<double>();
    params.monocyte_il_12_signal_virus_EC50 = (*cfg)["monocyte_il_12_signal_virus_EC50"].as<double>();
    params.monocyte_il_12_signal_virus_hill = (*cfg)["monocyte_il_12_signal_virus_hill"].as<double>();
    params.monocyte_il_12_downreg_il_10_EC50 = (*cfg)["monocyte_il_12_downreg_il_10_EC50"].as<double>();
    params.monocyte_il_12_downreg_il_10_hill = (*cfg)["monocyte_il_12_downreg_il_10_hill"].as<double>();

    /* G-CSF */
    params.monocyte_resting_basal_secreted_g_csf = (*cfg)["monocyte_resting_basal_secreted_g_csf"].as<double>();
    params.monocyte_resting_secreted_g_csf = (*cfg)["monocyte_resting_secreted_g_csf"].as<double>();
    params.monocyte_resting_uptake_g_csf = (*cfg)["monocyte_resting_uptake_g_csf"].as<double>();

    params.monocyte_active_basal_secreted_g_csf = (*cfg)["monocyte_active_basal_secreted_g_csf"].as<double>();
    params.monocyte_active_secreted_g_csf = (*cfg)["monocyte_active_secreted_g_csf"].as<double>();
    params.monocyte_active_uptake_g_csf = (*cfg)["monocyte_active_uptake_g_csf"].as<double>();

    params.monocyte_g_csf_EC50 = (*cfg)["monocyte_g_csf_EC50"].as<double>();
    params.monocyte_g_csf_hill = (*cfg)["monocyte_g_csf_hill"].as<double>();
    params.monocyte_g_csf_signal_virus_EC50 = (*cfg)["monocyte_g_csf_signal_virus_EC50"].as<double>();
    params.monocyte_g_csf_signal_virus_hill = (*cfg)["monocyte_g_csf_signal_virus_hill"].as<double>();
    params.monocyte_g_csf_downreg_il_10_EC50 = (*cfg)["monocyte_g_csf_downreg_il_10_EC50"].as<double>();
    params.monocyte_g_csf_downreg_il_10_hill = (*cfg)["monocyte_g_csf_downreg_il_10_hill"].as<double>();

    /* M-CSF */
    params.monocyte_resting_basal_secreted_m_csf = (*cfg)["monocyte_resting_basal_secreted_m_csf"].as<double>();
    params.monocyte_resting_secreted_m_csf = (*cfg)["monocyte_resting_secreted_m_csf"].as<double>();
    params.monocyte_resting_uptake_m_csf = (*cfg)["monocyte_resting_uptake_m_csf"].as<double>();

    params.monocyte_active_basal_secreted_m_csf = (*cfg)["monocyte_active_basal_secreted_m_csf"].as<double>();
    params.monocyte_active_secreted_m_csf = (*cfg)["monocyte_active_secreted_m_csf"].as<double>();
    params.monocyte_active_uptake_m_csf = (*cfg)["monocyte_active_uptake_m_csf"].as<double>();

    params.monocyte_m_csf_EC50 = (*cfg)["monocyte_m_csf_EC50"].as<double>();
    params.monocyte_m_csf_hill = (*cfg)["monocyte_m_csf_hill"].as<double>();
    params.monocyte_m_csf_signal_virus_EC50 = (*cfg)["monocyte_m_csf_signal_virus_EC50"].as<double>();
    params.monocyte_m_csf_signal_virus_hill = (*cfg)["monocyte_m_csf_signal_virus_hill"].as<double>();
    params.monocyte_m_csf_downreg_il_10_EC50 = (*cfg)["monocyte_m_csf_downreg_il_10_EC50"].as<double>();
    params.monocyte_m_csf_downreg_il_10_hill = (*cfg)["monocyte_m_csf_downreg_il_10_hill"].as<double>();

    /* infected epithelial chemokine */
    params.monocyte_resting_uptake_epi_inf_ck = (*cfg)["monocyte_resting_uptake_epi_inf_ck"].as<double>();
    params.monocyte_active_uptake_epi_inf_ck = (*cfg)["monocyte_active_uptake_epi_inf_ck"].as<double>();

    /* apoptosis chemokine */
    params.monocyte_secreted_apop_ck = (*cfg)["monocyte_secreted_apop_ck"].as<double>();

    /* phagocytosis chemokine */
    params.monocyte_resting_uptake_phago_ck = (*cfg)["monocyte_resting_uptake_phago_ck"].as<double>();
    params.monocyte_active_uptake_phago_ck = (*cfg)["monocyte_active_uptake_phago_ck"].as<double>();
    params.monocyte_apoptotic_secreted_phago_ck = (*cfg)["monocyte_apoptotic_secreted_phago_ck"].as<double>();

    /* infectious virus */
    params.monocyte_resting_uptake_virus_inf = (*cfg)["monocyte_resting_uptake_virus_inf"].as<double>();
    params.monocyte_resting_decay_virus_inf = (*cfg)["monocyte_resting_decay_virus_inf"].as<double>();

    params.monocyte_active_uptake_virus_inf = (*cfg)["monocyte_active_uptake_virus_inf"].as<double>();
    params.monocyte_active_decay_virus_inf = (*cfg)["monocyte_active_decay_virus_inf"].as<double>();

    params.monocyte_virus_inf_upreg_extracellular_virus_EC50 = (*cfg)["monocyte_virus_inf_upreg_extracellular_virus_EC50"].as<double>();
    params.monocyte_virus_inf_upreg_extracellular_virus_hill = (*cfg)["monocyte_virus_inf_upreg_extracellular_virus_hill"].as<double>();
    params.monocyte_virus_inf_downreg_internalised_virus_EC50 = (*cfg)["monocyte_virus_inf_downreg_internalised_virus_EC50"].as<double>();
    params.monocyte_virus_inf_downreg_internalised_virus_hill = (*cfg)["monocyte_virus_inf_downreg_internalised_virus_hill"].as<double>();

    /* non-infectious virus - no global memory; not used */

    /* movement */
    params.monocyte_resting_movement_rate = (*cfg)["monocyte_resting_movement_rate"].as<double>();
    params.monocyte_resting_movement_bias = (*cfg)["monocyte_resting_movement_bias"].as<double>();

    params.monocyte_active_movement_rate = (*cfg)["monocyte_active_movement_rate"].as<double>();
    params.monocyte_active_movement_bias = (*cfg)["monocyte_active_movement_bias"].as<double>();

    params.monocyte_movement_upreg_virus_EC50 = (*cfg)["monocyte_movement_upreg_virus_EC50"].as<double>();
    params.monocyte_movement_upreg_virus_hill = (*cfg)["monocyte_movement_upreg_virus_hill"].as<double>();

    params.monocyte_movement_weight_epi_inf_ck = (*cfg)["monocyte_movement_weight_epi_inf_ck"].as<double>();
    params.monocyte_movement_weight_phago_ck = (*cfg)["monocyte_movement_weight_phago_ck"].as<double>();

    /* recruitment */
    params.monocyte_recruitment_weight_il_6 = (*cfg)["monocyte_recruitment_weight_il_6"].as<double>();
    params.monocyte_recruitment_weight_il_10 = (*cfg)["monocyte_recruitment_weight_il_10"].as<double>();
    params.monocyte_recruitment_weight_m_csf = (*cfg)["monocyte_recruitment_weight_m_csf"].as<double>();
    params.monocyte_recruitment_weight_epi_inf_ck = (*cfg)["monocyte_recruitment_weight_epi_inf_ck"].as<double>();
    params.monocyte_recruitment_weight_phago_ck = (*cfg)["monocyte_recruitment_weight_phago_ck"].as<double>();

    params.monocyte_recruitment_EC50 = (*cfg)["monocyte_recruitment_EC50"].as<double>();
    params.monocyte_recruitment_hill = (*cfg)["monocyte_recruitment_hill"].as<double>();

    /* activation and deactivation */
    params.monocyte_activate_weight_ifn_1 = (*cfg)["monocyte_activate_weight_ifn_1"].as<double>();
    params.monocyte_activate_weight_il_6 = (*cfg)["monocyte_activate_weight_il_6"].as<double>();
    params.monocyte_activate_weight_il_10 = (*cfg)["monocyte_activate_weight_il_10"].as<double>();
    params.monocyte_activate_weight_m_csf = (*cfg)["monocyte_activate_weight_m_csf"].as<double>();
    params.monocyte_activate_weight_epi_inf_ck = (*cfg)["monocyte_activate_weight_epi_inf_ck"].as<double>();
    params.monocyte_activate_weight_phago_ck = (*cfg)["monocyte_activate_weight_phago_ck"].as<double>();
    params.monocyte_activate_weight_virus_inf = (*cfg)["monocyte_activate_weight_virus_inf"].as<double>();

    params.monocyte_activate_EC50 = (*cfg)["monocyte_activate_EC50"].as<double>();
    params.monocyte_activate_hill = (*cfg)["monocyte_activate_hill"].as<double>();

    /* bursting */
    params.monocyte_resting_burst_threshold = (*cfg)["monocyte_resting_burst_threshold"].as<double>();
    params.monocyte_active_burst_threshold = (*cfg)["monocyte_active_burst_threshold"].as<double>();

    /* apoptosis */
    params.monocyte_apoptosis_timescale = (*cfg)["monocyte_apoptosis_timescale"].as<double>();

    /*==========*/
    /* not used */
    params.monocyte_apoptosis_max_rate = (*cfg)["monocyte_apoptosis_max_rate"].as<double>();
    params.monocyte_apoptosis_half_max = (*cfg)["monocyte_apoptosis_half_max"].as<double>();
    params.monocyte_apoptosis_hill = (*cfg)["monocyte_apoptosis_hill"].as<double>();
    /*==========*/

    params.monocyte_apoptosis_by_ifn_max_rate = (*cfg)["monocyte_apoptosis_by_ifn_max_rate"].as<double>();
    params.monocyte_apoptosis_by_ifn_half_max = (*cfg)["monocyte_apoptosis_by_ifn_half_max"].as<double>();
    params.monocyte_apoptosis_by_ifn_hill = (*cfg)["monocyte_apoptosis_by_ifn_hill"].as<double>();

    /* phagocytosis */
    params.monocyte_phagocytosis_timescale = (*cfg)["monocyte_phagocytosis_timescale"].as<double>();
    params.monocyte_active_phagocytose_apoptotic_epithelial_probability = (*cfg)["monocyte_active_phagocytose_apoptotic_epithelial_probability"].as<double>();
    params.monocyte_active_phagocytose_infected_epithelial_probability = (*cfg)["monocyte_active_phagocytose_infected_epithelial_probability"].as<double>();
    params.monocyte_active_phagocytose_immune_probability = (*cfg)["monocyte_active_phagocytose_immune_probability"].as<double>();

    /* conversion */
    params.monocyte_active_to_mac_probability = (*cfg)["monocyte_active_to_mac_probability"].as<double>();

    params.monocyte_blocking_function = false;

  } catch (std::exception &e) {
    CAM_ERROR(e.what(), CAM_ERROR_IO);
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Monocyte constructor
 * =============================================================================
 */
Monocyte::Monocyte (const int x, const int y, const MonocyteType tp, const MonocyteState st,
                    struct cam_rd_sys* rds, const bool use_output, int* status)
  : Agent{ x, y, AgentType::monocyte, rds }
{
  /* agent is constructed first */

  /* set monocyte type */
  set_type(tp);

  MonocyteState state;
  if (st == MonocyteState::null) {
    /* calculate random number between 0 and 1 */
    double r = distribution_0_1(rng);
    if (r < params.monocyte_active_probability) {
      /* monocyte state is active */
      state = MonocyteState::active;
    } else {
      /* monocyte state is resting */
      state = MonocyteState::resting;
    }
  } else {
    state = st;
  }

  total_virus_uptaken = 0.0;

  /* set monocyte state, age, lifespan, secretions and uptake */
  set_state(state);

  /* set the area */
  double mac_area = CAM_PI * params.monocyte_radius * params.monocyte_radius;
  set_area(mac_area);

  /* update the reaction diffusion secretions & uptakes */
  inc_rds_values();

  /* resident on cell */
  resident_on_cell = nullptr;

  /* apoptosis */
  is_being_apoptosed = false;
  apoptosis_start_time = 0.0;
  apoptosis_timescale = params.monocyte_apoptosis_timescale;
  num_cells_apoptosing = 0;
  is_apoptosing_epithelial = false;
  apoptosing_agent = nullptr;

  /* phagocytosis */
  phagocytosis_timescale = params.monocyte_phagocytosis_timescale;
  is_phagocytosing_epithelial = false;
  phagocytosing_agent = nullptr;
  phagocytosed_by = nullptr;

  std::string outf;
  if (use_output) {
    outf = output_folder + std::string("/path_monocyte.txt");
    output = fopen(outf.c_str(), "w");

    outf = output_folder + std::string("/monocyte_total_virus_uptaken.txt");
    uptake_output = fopen(outf.c_str(), "w");

    outf = output_folder + std::string("/monocyte_lifespan.txt");
    lifespan_output = fopen(outf.c_str(), "w");

    outf = output_folder + std::string("/monocyte_movement_rate.txt");
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
 * Monocyte destructor
 * =============================================================================
 */
Monocyte::~Monocyte ()
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
 * set monocyte state, reset age, lifespan, secretions and uptakes
 * =============================================================================
 */
void
Monocyte::set_state (const MonocyteState st)
{
  /* set state */
  mono_state = st;

  /* reset age */
  set_age(0.0);

  /* set the lifespan, depending on type and state */
  double lspan = 0.0;
  switch (mono_state) {
    case MonocyteState::resting:
      lspan = params.monocyte_resting_lifespan;
      break;
    case MonocyteState::active:
      lspan = params.monocyte_active_lifespan;
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
Monocyte::set_secretion (const MonocyteState st)
{
  switch (st) {
    case MonocyteState::resting:
      basal_secreted_ifn_1 = params.monocyte_resting_basal_secreted_ifn_1;
      secreted_ifn_1 = params.monocyte_resting_secreted_ifn_1;

      basal_secreted_il_6 = params.monocyte_resting_basal_secreted_il_6;
      secreted_il_6 = params.monocyte_resting_secreted_il_6;

      basal_secreted_il_10 = params.monocyte_resting_basal_secreted_il_10;
      secreted_il_10 = params.monocyte_resting_secreted_il_10;

      basal_secreted_il_12 = params.monocyte_resting_basal_secreted_il_12;
      secreted_il_12 = params.monocyte_resting_secreted_il_12;

      basal_secreted_g_csf = params.monocyte_resting_basal_secreted_g_csf;
      secreted_g_csf = params.monocyte_resting_secreted_g_csf;

      basal_secreted_m_csf = params.monocyte_resting_basal_secreted_m_csf;
      secreted_m_csf = params.monocyte_resting_secreted_m_csf;

      basal_secreted_epi_inf_ck = 0.0;
      secreted_epi_inf_ck = 0.0;

      secreted_apop_ck = params.monocyte_secreted_apop_ck;
      secreted_phago_ck = 0.0;

      break;
    case MonocyteState::active:
      basal_secreted_ifn_1 = params.monocyte_active_basal_secreted_ifn_1;
      secreted_ifn_1 = params.monocyte_active_secreted_ifn_1;

      basal_secreted_il_6 = params.monocyte_active_basal_secreted_il_6;
      secreted_il_6 = params.monocyte_active_secreted_il_6;

      basal_secreted_il_10 = params.monocyte_active_basal_secreted_il_10;
      secreted_il_10 = params.monocyte_active_secreted_il_10;

      basal_secreted_il_12 = params.monocyte_active_basal_secreted_il_12;
      secreted_il_12 = params.monocyte_active_secreted_il_12;

      basal_secreted_g_csf = params.monocyte_active_basal_secreted_g_csf;
      secreted_g_csf = params.monocyte_active_secreted_g_csf;

      basal_secreted_m_csf = params.monocyte_active_basal_secreted_m_csf;
      secreted_m_csf = params.monocyte_active_secreted_m_csf;

      basal_secreted_epi_inf_ck = 0.0;
      secreted_epi_inf_ck = 0.0;

      secreted_apop_ck = params.monocyte_secreted_apop_ck;
      secreted_phago_ck = 0.0;

      break;
    case MonocyteState::apoptotic:
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
      secreted_phago_ck = params.monocyte_apoptotic_secreted_phago_ck;

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
Monocyte::set_uptake (const MonocyteState st)
{
  switch (st) {
    case MonocyteState::resting:
      uptake_ifn_1 = params.monocyte_resting_uptake_ifn_1;

      uptake_il_6 = params.monocyte_resting_uptake_il_6;
      uptake_il_10 = params.monocyte_resting_uptake_il_10;
      uptake_il_12 = params.monocyte_resting_uptake_il_12;

      uptake_g_csf = params.monocyte_resting_uptake_g_csf;
      uptake_m_csf = params.monocyte_resting_uptake_m_csf;

      uptake_epi_inf_ck = params.monocyte_resting_uptake_epi_inf_ck;
      uptake_apop_ck = 0.0;
      uptake_phago_ck = params.monocyte_resting_uptake_phago_ck;

      uptake_virus_inf = params.monocyte_resting_uptake_virus_inf;
      decay_virus_inf = params.monocyte_resting_decay_virus_inf;

      break;
    case MonocyteState::active:
      uptake_ifn_1 = params.monocyte_active_uptake_ifn_1;

      uptake_il_6 = params.monocyte_active_uptake_il_6;
      uptake_il_10 = params.monocyte_active_uptake_il_10;
      uptake_il_12 = params.monocyte_active_uptake_il_12;

      uptake_g_csf = params.monocyte_active_uptake_g_csf;
      uptake_m_csf = params.monocyte_active_uptake_m_csf;

      uptake_epi_inf_ck = params.monocyte_active_uptake_epi_inf_ck;
      uptake_apop_ck = 0.0;
      uptake_phago_ck = params.monocyte_active_uptake_phago_ck;

      uptake_virus_inf = params.monocyte_active_uptake_virus_inf;
      decay_virus_inf = params.monocyte_active_decay_virus_inf;

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
Monocyte::set_movement (const MonocyteState st)
{
  switch (st) {
    case MonocyteState::resting:
      movement_rate = params.monocyte_resting_movement_rate;
      movement_bias = params.monocyte_resting_movement_bias;
      break;
    case MonocyteState::active:
      movement_rate = params.monocyte_active_movement_rate;
      movement_bias = params.monocyte_active_movement_bias;
      break;
    default:
      movement_rate = 0.0;
      movement_bias = 0.0;
      break;
  }
}


void
Monocyte::change_state (const MonocyteState st)
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
Monocyte::update_lifespan (const MonocyteState st)
{
  (void) st;

  return;

  double inhib = cam_util_hill_function(
    params.monocyte_lifespan_downreg_virus_EC50, total_virus_uptaken, params.monocyte_lifespan_downreg_virus_hill
  );

  double lspan = get_lifespan(true)*inhib; // get_lifespan(true) returns max lifespan
  set_lifespan(lspan, false);
}


void
Monocyte::update_movement (const MonocyteState st)
{
  /* reset the movement rate */
  set_movement(st);

  return;

  double conc = cam_util_hill_function(
    total_virus_uptaken, params.monocyte_movement_upreg_virus_EC50, params.monocyte_movement_upreg_virus_hill
  );

  movement_rate *= (1.0 + conc); // NEEDSWORK: => max at 2*movement_rate
}


void
Monocyte::regulate_secretion ()
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
  conc = total_virus_uptaken;
  signal = cam_util_hill_function(
    conc, params.monocyte_ifn_1_signal_virus_EC50, params.monocyte_ifn_1_signal_virus_hill
  );

  conc = il_10;
  downreg = cam_util_hill_function(
    params.monocyte_ifn_1_downreg_il_10_EC50, conc, params.monocyte_ifn_1_downreg_il_10_hill
  );

  conc = ifn_1;
  upreg = cam_util_hill_function(
    conc, params.monocyte_ifn_1_EC50, params.monocyte_ifn_1_hill
  );

  secreted_ifn_1 = (
    ((basal_secreted_ifn_1 / (1.0 + conc*conc)) + secreted_ifn_1*upreg*downreg)*signal
  );


  /* IL-6 */
  conc = total_virus_uptaken;
  signal = cam_util_hill_function(
    conc, params.monocyte_il_6_signal_virus_EC50, params.monocyte_il_6_signal_virus_hill
  );

  conc = il_10;
  downreg = cam_util_hill_function(
    params.monocyte_il_6_downreg_il_10_EC50, conc, params.monocyte_il_6_downreg_il_10_hill
  );

  conc = il_6;
  upreg = cam_util_hill_function(
    conc, params.monocyte_il_6_EC50, params.monocyte_il_6_hill
  );

  secreted_il_6 = (
    ((basal_secreted_il_6 / (1.0 + conc*conc)) + secreted_il_6*upreg*downreg)*signal
  );


  /* IL-12 */
  conc = total_virus_uptaken;
  signal = cam_util_hill_function(
    conc, params.monocyte_il_12_signal_virus_EC50, params.monocyte_il_12_signal_virus_hill
  );

  conc = il_10;
  downreg = cam_util_hill_function(
    params.monocyte_il_12_downreg_il_10_EC50, conc, params.monocyte_il_12_downreg_il_10_hill
  );

  conc = il_12;
  upreg = cam_util_hill_function(
    conc, params.monocyte_il_12_EC50, params.monocyte_il_12_hill
  );

  secreted_il_12 = (
    ((basal_secreted_il_12 / (1.0 + conc*conc)) + secreted_il_12*upreg*downreg)*signal
  );


  /* IL-10 */
  conc = il_6;
  signal = cam_util_hill_function(
    conc, params.monocyte_il_10_signal_il_6_EC50, params.monocyte_il_10_signal_il_6_hill
  );
  conc = il_12;
  signal += cam_util_hill_function(
    conc, params.monocyte_il_10_signal_il_12_EC50, params.monocyte_il_10_signal_il_12_hill
  );
  conc = g_csf;
  signal += cam_util_hill_function(
    conc, params.monocyte_il_10_signal_g_csf_EC50, params.monocyte_il_10_signal_g_csf_hill
  );
  conc = m_csf;
  signal += cam_util_hill_function(
    conc, params.monocyte_il_10_signal_m_csf_EC50, params.monocyte_il_10_signal_m_csf_hill
  );

  conc = il_10;
  upreg = cam_util_hill_function(
    conc, params.monocyte_il_10_EC50, params.monocyte_il_10_hill
  );

  secreted_il_10 = (
    ((basal_secreted_il_10 / (1.0 + conc*conc)) + secreted_il_10*upreg)*signal
  );


  /* G-CSF */
  conc = total_virus_uptaken;
  signal = cam_util_hill_function(
    conc, params.monocyte_g_csf_signal_virus_EC50, params.monocyte_g_csf_signal_virus_hill
  );

  conc = il_10;
  downreg = cam_util_hill_function(
    params.monocyte_g_csf_downreg_il_10_EC50, conc, params.monocyte_g_csf_downreg_il_10_hill
  );

  conc = g_csf;
  upreg = cam_util_hill_function(
    conc, params.monocyte_g_csf_EC50, params.monocyte_g_csf_hill
  );

  secreted_g_csf = (
    ((basal_secreted_g_csf / (1.0 + conc*conc)) + secreted_g_csf*upreg*downreg)*signal
  );


  /* M-CSF */
  conc = total_virus_uptaken;
  signal = cam_util_hill_function(
    conc, params.monocyte_m_csf_signal_virus_EC50, params.monocyte_m_csf_signal_virus_hill
  );

  conc = il_10;
  downreg = cam_util_hill_function(
    params.monocyte_m_csf_downreg_il_10_EC50, conc, params.monocyte_m_csf_downreg_il_10_hill
  );

  conc = m_csf;
  upreg = cam_util_hill_function(
    conc, params.monocyte_m_csf_EC50, params.monocyte_m_csf_hill
  );

  secreted_m_csf = (
    ((basal_secreted_m_csf / (1.0 + conc*conc)) + secreted_m_csf*upreg*downreg)*signal
  );


  /* infected epithelial chemokine - not used */
}


void
Monocyte::regulate_uptake ()
{
  std::pair<int, int> xy = get_coordinates();

  double epi_area = ( (resident_on_cell == nullptr) ? 0.0 : resident_on_cell->get_max_area());
  double n_v = rd_sys->virus_inf->get_conc(xy.first, xy.second) * epi_area;

  double gfnc = 0.0;
  if (params.monocyte_blocking_function) {
    gfnc = 1.0;
  } else if (n_v >= 1.0) {
    params.monocyte_blocking_function = true;
    gfnc = 1.0;
  }

  uptake_virus_inf *= n_v * gfnc * cam_util_hill_function(
    params.monocyte_virus_inf_downreg_internalised_virus_EC50,
    total_virus_uptaken,
    params.monocyte_virus_inf_downreg_internalised_virus_hill
  );
}


void
Monocyte::update_secretions ()
{
  /* update the reaction diffusion secretions and uptakes */
  dec_rds_values();

  /* update the lifespan */
  update_lifespan(mono_state);

  /* reset the secretions */
  set_secretion(mono_state);

  /* reset the uptakes */
  set_uptake(mono_state);

  /* update the movement */
  update_movement(mono_state);

  /* update the reaction diffusion secretions and uptakes */
  inc_rds_values();
}


/*
 * =============================================================================
 * Compute the monocyte action (activation, deactivation, phagocytosis)
 * =============================================================================
 */
void
Monocyte::action (const double time)
{
  /* random number */
  double r;

  /* coords */
  std::pair<int, int> xy = get_coordinates();

  if (mono_state == MonocyteState::apoptotic || is_being_apoptosed || is_being_phagocytosed || mono_state == MonocyteState::phagocytosed) {
    return;
  }

  /* update secretions */
  update_secretions();

  double ifn_uptake = rd_sys->ifn_1->get_conc(xy.first, xy.second);
  double trigger_apop = params.monocyte_apoptosis_by_ifn_max_rate * cam_util_hill_function(ifn_uptake, params.monocyte_apoptosis_by_ifn_half_max, params.monocyte_apoptosis_by_ifn_hill);

  r = distribution_0_1(rng);
  if (r < trigger_apop) {
    increase_age(10.0*get_lifespan(true));

    return;
  }

  total_virus_uptaken += timestep * (uptake_virus_inf - decay_virus_inf * total_virus_uptaken);

  /* get the aggregate of the cytokines */
  double cyt_tmp = (
    params.monocyte_activate_weight_virus_inf * total_virus_uptaken +
    params.monocyte_activate_weight_ifn_1 * rd_sys->ifn_1->get_conc(xy.first, xy.second) +
    params.monocyte_activate_weight_il_6 * rd_sys->il_6->get_conc(xy.first, xy.second) +
    params.monocyte_activate_weight_m_csf * rd_sys->m_csf->get_conc(xy.first, xy.second) +
    params.monocyte_activate_weight_epi_inf_ck * rd_sys->epi_inf_ck->get_conc(xy.first, xy.second) +
    params.monocyte_activate_weight_phago_ck * rd_sys->phago_ck->get_conc(xy.first, xy.second) -
    params.monocyte_activate_weight_il_10 * rd_sys->il_10->get_conc(xy.first, xy.second)
  );

  double cyt = cam_util_hill_function(
    cyt_tmp, params.monocyte_activate_EC50, params.monocyte_activate_hill
  );

  /* get the area */
  const double mono_area = get_area();

  /* resting actions */
  if (mono_state == MonocyteState::resting) {
    /* if total_virus_uptaken is beyond the threshold then burst */
    if (total_virus_uptaken*mono_area > params.monocyte_resting_burst_threshold) {
      /* flag for bursting - handled in environment */
      flag_bursting = true;
      return;
    }

    /* activation */
    r = distribution_0_1(rng);
    if (r < cyt && !is_phagocytosing_epithelial && phagocytosing_agent == nullptr && !is_being_apoptosed && !is_being_phagocytosed) {
      change_state(MonocyteState::active);

      /* update the secretions */
      update_secretions();

      return;
    }
  }


  /* active actions */
  if (mono_state == MonocyteState::active) {
    /* if total_virus_uptaken is beyond the threshold then burst */
    if (total_virus_uptaken*mono_area > params.monocyte_active_burst_threshold) {
      /* flag for bursting - handled in environment */
      flag_bursting = true;
      return;
    }

    /* conversion promoting cytokines */
    double cyt_promote = (
      rd_sys->m_csf->get_scaled(xy.first, xy.second) + rd_sys->il_6->get_scaled(xy.first, xy.second)
    ) / 2.0;

    /* flag for conversion */
    r = distribution_0_1(rng);
    if (r < params.monocyte_active_to_mac_probability*(1.0 + cyt_promote)) {
      this->flag_conversion = true; // protected element of Agent class

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
    if (r > cyt && !is_phagocytosing_epithelial && phagocytosing_agent == nullptr && !is_being_apoptosed && !is_being_phagocytosed) {
      change_state(MonocyteState::resting);

      /* update the secretions */
      update_secretions();

      return;
    }
  }
}


int
Monocyte::phagocytose_epithelial (const double time)
{
  if (phagocytosing_agent != nullptr) {
    /* this monocyte cell is already phagocytosing another immune cell */
    return -1;
  }


  /* phagocytosis of epithelial cell */
  if (is_phagocytosing_epithelial && !resident_on_cell->is_being_phagocytosed) {
    /* phagocytosis is complete */

    is_phagocytosing_epithelial = false;

    return 0;

  } else if (!is_phagocytosing_epithelial && !resident_on_cell->is_being_phagocytosed) {
    /* monocyte starting phagocytosis */

    if (resident_on_cell->check_state(EpithelialState::apoptotic)) {
      double r = distribution_0_1(rng);
      if (r < params.monocyte_active_phagocytose_apoptotic_epithelial_probability) {
        /* start phagocytosis */

        double tscale = phagocytosis_timescale;

        resident_on_cell->start_phagocytosis(time);
        resident_on_cell->set_phagocytosis_timescale(tscale);
        resident_on_cell->switch_flag_phagocytosis();
        is_phagocytosing_epithelial = true;

        return 0;
      }
    } else if (resident_on_cell->check_state(EpithelialState::infected)) {
      double r = distribution_0_1(rng);
      if (r < params.monocyte_active_phagocytose_infected_epithelial_probability) {
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
Monocyte::phagocytose_immune (const double time)
{
  if (is_phagocytosing_epithelial) {
    /* this monocyte is phagocytosing an epithelial cell then return */
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
      if (r < (params.monocyte_active_phagocytose_immune_probability)) {
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
Monocyte::migrate (const double dt, const double tm, std::function<int(Agent*, std::vector<std::pair<int,int>>)> fndsp, std::function<bool(Agent*, std::pair<int,int>)> chksp)
{
  Epithelial* epi = resident_on_cell;
  const int sz = (int) epi->immediate_moore.size();

  /* return index */
  int ind = sz;

  /* movement rate */
  if (is_phagocytosing_epithelial || phagocytosing_agent != nullptr || is_being_phagocytosed ||
      is_being_apoptosed || (mono_state == MonocyteState::apoptotic)) {
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
        params.monocyte_movement_weight_epi_inf_ck * rd_sys->epi_inf_ck->get_conc(xy.first, xy.second) +
        params.monocyte_movement_weight_phago_ck * rd_sys->phago_ck->get_conc(xy.first, xy.second)
      );

      total_cyt += pow(nbr_cyt, movement_bias);
    }

    if (total_cyt < CAM_TOL || mono_state == MonocyteState::resting) {
      /* RANDOM MIGRATION: choose an immediate neighbour */

      /* find space using function provided as arg */
      ind = fndsp(this, epi->immediate_moore);
    } else {
      /* CHEMOTACTIC MIGRATION: (biased random walk) choose neighbour based on cytokine levels */

      double r = distribution_0_1(rng)*total_cyt;
      double running_total = 0.0;
      for (int j = 0; j < sz; j++) {
        /* coordinates of neighbour */
        xy = epi->immediate_moore[j];

        /* cytokine/chemokine concentration used for migration */
        nbr_cyt = (
          params.monocyte_movement_weight_epi_inf_ck * rd_sys->epi_inf_ck->get_conc(xy.first, xy.second) +
          params.monocyte_movement_weight_phago_ck * rd_sys->phago_ck->get_conc(xy.first, xy.second)
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
Monocyte::recruit (std::vector<std::pair<int,int>> nbrs)
{
  if (params.monocyte_recruitment_hill < CAM_TOL) {
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
    params.monocyte_recruitment_weight_m_csf * nbr_m_csf +
    params.monocyte_recruitment_weight_il_6 * nbr_il_6 +
    params.monocyte_recruitment_weight_epi_inf_ck * nbr_epi_inf_ck +
    params.monocyte_recruitment_weight_phago_ck * nbr_phago_ck -
    params.monocyte_recruitment_weight_il_10 * nbr_il_10
  ) / ((double) nbrs.size());

  double cyt = cam_util_hill_function(
    total_cyt, params.monocyte_recruitment_EC50, params.monocyte_recruitment_hill
  );

  double r = distribution_0_1(rng);
  if (r < cyt) {
    return true;
  }

  return false;
}


/*
 * =============================================================================
 * Start the apoptosis in the monocyte
 * =============================================================================
 */
void
Monocyte::start_apoptosis (const double time)
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
 * End the apoptosis in the monocyte
 * =============================================================================
 */
void
Monocyte::end_apoptosis ()
{
  std::pair<int, int> xy = get_coordinates();

  is_being_apoptosed = false;
  apoptosis_start_time = 0.0;
  num_cells_apoptosing = 0;

  change_state(MonocyteState::apoptotic);

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
Monocyte::set_apoptosis_timescale (const double tms)
{
  double timescale = tms / ((double) num_cells_apoptosing);

  /* change time scale of apoptosis due to NK action */
  if (timescale < apoptosis_timescale) {
    apoptosis_timescale = timescale;
  }
}


/*
 * =============================================================================
 * Start the phagocytosis in the monocyte
 * =============================================================================
 */
void
Monocyte::start_phagocytosis (const double time)
{
  is_being_phagocytosed = true;
  phagocytosis_start_time = time;
}


/*
 * =============================================================================
 * End the phagocytosis in the monocyte
 * =============================================================================
 */
void
Monocyte::end_phagocytosis ()
{
  is_being_phagocytosed = false;
  phagocytosis_start_time = 0.0;

  phagocytosed_by->phagocytosing_agent = nullptr;
  phagocytosed_by = nullptr;
  phagocytosing_agent = nullptr;

  change_state(MonocyteState::phagocytosed);

  /* update the secretions */
  update_secretions();
}


/*
 * =============================================================================
 * Set the phagocytosis timescale
 * =============================================================================
 */
void
Monocyte::set_phagocytosis_timescale (const double tms)
{
  (void) tms; // phagocytosis is set by the target cell only

  return;
  
  /* change time scale of phagocytosis */
  if (tms < phagocytosis_timescale) {
    phagocytosis_timescale = tms;
  }
}
