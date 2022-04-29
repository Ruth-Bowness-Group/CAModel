#include <cstring>

#include "cam_nk.h"
#include "nk.h"
#include "../../cam_error.h"
#include "../../cam_random.h"
#include "../../cam_util.h"
#include "../../cam_parameters.h"


/* struct containing nk specific parameters */
static nk_params params;


/*
 * =============================================================================
 * Set the nk specific parameters
 * =============================================================================
 */
int
cam_nk_set_params (YAML::Node* cfg)
{
  std::string ren;

  try {
    /* radius */
    params.nk_radius = (*cfg)["nk_radius"].as<double>();

    /* probability of activating nk */
    params.nk_active_probability = (*cfg)["nk_active_probability"].as<double>();

    /* lifespan */
    params.nk_resting_lifespan = (*cfg)["nk_resting_lifespan"].as<double>();
    params.nk_active_lifespan = (*cfg)["nk_active_lifespan"].as<double>();

    params.nk_lifespan_downreg_virus_EC50 = (*cfg)["nk_lifespan_downreg_virus_EC50"].as<double>();
    params.nk_lifespan_downreg_virus_hill = (*cfg)["nk_lifespan_downreg_virus_hill"].as<double>();

    /*
     * cytokines
     */

    /* T1 IFN */
    params.nk_resting_basal_secreted_ifn_1 = (*cfg)["nk_resting_basal_secreted_ifn_1"].as<double>();
    params.nk_resting_secreted_ifn_1 = (*cfg)["nk_resting_secreted_ifn_1"].as<double>();
    params.nk_resting_uptake_ifn_1 = (*cfg)["nk_resting_uptake_ifn_1"].as<double>();

    params.nk_active_basal_secreted_ifn_1 = (*cfg)["nk_active_basal_secreted_ifn_1"].as<double>();
    params.nk_active_secreted_ifn_1 = (*cfg)["nk_active_secreted_ifn_1"].as<double>();
    params.nk_active_uptake_ifn_1 = (*cfg)["nk_active_uptake_ifn_1"].as<double>();

    params.nk_ifn_1_EC50 = (*cfg)["nk_ifn_1_EC50"].as<double>();
    params.nk_ifn_1_hill = (*cfg)["nk_ifn_1_hill"].as<double>();
    params.nk_ifn_1_signal_virus_EC50 = (*cfg)["nk_ifn_1_signal_virus_EC50"].as<double>();
    params.nk_ifn_1_signal_virus_hill = (*cfg)["nk_ifn_1_signal_virus_hill"].as<double>();
    params.nk_ifn_1_downreg_il_10_EC50 = (*cfg)["nk_ifn_1_downreg_il_10_EC50"].as<double>();
    params.nk_ifn_1_downreg_il_10_hill = (*cfg)["nk_ifn_1_downreg_il_10_hill"].as<double>();

    /* IL-6 */
    params.nk_resting_basal_secreted_il_6 = (*cfg)["nk_resting_basal_secreted_il_6"].as<double>();
    params.nk_resting_secreted_il_6 = (*cfg)["nk_resting_secreted_il_6"].as<double>();
    params.nk_resting_uptake_il_6 = (*cfg)["nk_resting_uptake_il_6"].as<double>();

    params.nk_active_basal_secreted_il_6 = (*cfg)["nk_active_basal_secreted_il_6"].as<double>();
    params.nk_active_secreted_il_6 = (*cfg)["nk_active_secreted_il_6"].as<double>();
    params.nk_active_uptake_il_6 = (*cfg)["nk_active_uptake_il_6"].as<double>();

    params.nk_il_6_EC50 = (*cfg)["nk_il_6_EC50"].as<double>();
    params.nk_il_6_hill = (*cfg)["nk_il_6_hill"].as<double>();
    params.nk_il_6_signal_virus_EC50 = (*cfg)["nk_il_6_signal_virus_EC50"].as<double>();
    params.nk_il_6_signal_virus_hill = (*cfg)["nk_il_6_signal_virus_hill"].as<double>();
    params.nk_il_6_downreg_il_10_EC50 = (*cfg)["nk_il_6_downreg_il_10_EC50"].as<double>();
    params.nk_il_6_downreg_il_10_hill = (*cfg)["nk_il_6_downreg_il_10_hill"].as<double>();

    /* IL-10 */
    params.nk_resting_basal_secreted_il_10 = (*cfg)["nk_resting_basal_secreted_il_10"].as<double>();
    params.nk_resting_secreted_il_10 = (*cfg)["nk_resting_secreted_il_10"].as<double>();
    params.nk_resting_uptake_il_10 = (*cfg)["nk_resting_uptake_il_10"].as<double>();

    params.nk_active_basal_secreted_il_10 = (*cfg)["nk_active_basal_secreted_il_10"].as<double>();
    params.nk_active_secreted_il_10 = (*cfg)["nk_active_secreted_il_10"].as<double>();
    params.nk_active_uptake_il_10 = (*cfg)["nk_active_uptake_il_10"].as<double>();

    params.nk_il_10_EC50 = (*cfg)["nk_il_10_EC50"].as<double>();
    params.nk_il_10_hill = (*cfg)["nk_il_10_hill"].as<double>();
    params.nk_il_10_signal_il_6_EC50 = (*cfg)["nk_il_10_signal_il_6_EC50"].as<double>();
    params.nk_il_10_signal_il_6_hill = (*cfg)["nk_il_10_signal_il_6_hill"].as<double>();
    params.nk_il_10_signal_il_12_EC50 = (*cfg)["nk_il_10_signal_il_12_EC50"].as<double>();
    params.nk_il_10_signal_il_12_hill = (*cfg)["nk_il_10_signal_il_12_hill"].as<double>();
    params.nk_il_10_signal_g_csf_EC50 = (*cfg)["nk_il_10_signal_g_csf_EC50"].as<double>();
    params.nk_il_10_signal_g_csf_hill = (*cfg)["nk_il_10_signal_g_csf_hill"].as<double>();
    params.nk_il_10_signal_m_csf_EC50 = (*cfg)["nk_il_10_signal_m_csf_EC50"].as<double>();
    params.nk_il_10_signal_m_csf_hill = (*cfg)["nk_il_10_signal_m_csf_hill"].as<double>();

    /* IL-12 */
    params.nk_resting_basal_secreted_il_12 = (*cfg)["nk_resting_basal_secreted_il_12"].as<double>();
    params.nk_resting_secreted_il_12 = (*cfg)["nk_resting_secreted_il_12"].as<double>();
    params.nk_resting_uptake_il_12 = (*cfg)["nk_resting_uptake_il_12"].as<double>();

    params.nk_active_basal_secreted_il_12 = (*cfg)["nk_active_basal_secreted_il_12"].as<double>();
    params.nk_active_secreted_il_12 = (*cfg)["nk_active_secreted_il_12"].as<double>();
    params.nk_active_uptake_il_12 = (*cfg)["nk_active_uptake_il_12"].as<double>();

    params.nk_il_12_EC50 = (*cfg)["nk_il_12_EC50"].as<double>();
    params.nk_il_12_hill = (*cfg)["nk_il_12_hill"].as<double>();
    params.nk_il_12_signal_virus_EC50 = (*cfg)["nk_il_12_signal_virus_EC50"].as<double>();
    params.nk_il_12_signal_virus_hill = (*cfg)["nk_il_12_signal_virus_hill"].as<double>();
    params.nk_il_12_downreg_il_10_EC50 = (*cfg)["nk_il_12_downreg_il_10_EC50"].as<double>();
    params.nk_il_12_downreg_il_10_hill = (*cfg)["nk_il_12_downreg_il_10_hill"].as<double>();

    /* G-CSF */
    params.nk_resting_basal_secreted_g_csf = (*cfg)["nk_resting_basal_secreted_g_csf"].as<double>();
    params.nk_resting_secreted_g_csf = (*cfg)["nk_resting_secreted_g_csf"].as<double>();
    params.nk_resting_uptake_g_csf = (*cfg)["nk_resting_uptake_g_csf"].as<double>();

    params.nk_active_basal_secreted_g_csf = (*cfg)["nk_active_basal_secreted_g_csf"].as<double>();
    params.nk_active_secreted_g_csf = (*cfg)["nk_active_secreted_g_csf"].as<double>();
    params.nk_active_uptake_g_csf = (*cfg)["nk_active_uptake_g_csf"].as<double>();

    params.nk_g_csf_EC50 = (*cfg)["nk_g_csf_EC50"].as<double>();
    params.nk_g_csf_hill = (*cfg)["nk_g_csf_hill"].as<double>();
    params.nk_g_csf_signal_virus_EC50 = (*cfg)["nk_g_csf_signal_virus_EC50"].as<double>();
    params.nk_g_csf_signal_virus_hill = (*cfg)["nk_g_csf_signal_virus_hill"].as<double>();
    params.nk_g_csf_downreg_il_10_EC50 = (*cfg)["nk_g_csf_downreg_il_10_EC50"].as<double>();
    params.nk_g_csf_downreg_il_10_hill = (*cfg)["nk_g_csf_downreg_il_10_hill"].as<double>();

    /* M-CSF */
    params.nk_resting_basal_secreted_m_csf = (*cfg)["nk_resting_basal_secreted_m_csf"].as<double>();
    params.nk_resting_secreted_m_csf = (*cfg)["nk_resting_secreted_m_csf"].as<double>();
    params.nk_resting_uptake_m_csf = (*cfg)["nk_resting_uptake_m_csf"].as<double>();

    params.nk_active_basal_secreted_m_csf = (*cfg)["nk_active_basal_secreted_m_csf"].as<double>();
    params.nk_active_secreted_m_csf = (*cfg)["nk_active_secreted_m_csf"].as<double>();
    params.nk_active_uptake_m_csf = (*cfg)["nk_active_uptake_m_csf"].as<double>();

    params.nk_m_csf_EC50 = (*cfg)["nk_m_csf_EC50"].as<double>();
    params.nk_m_csf_hill = (*cfg)["nk_m_csf_hill"].as<double>();
    params.nk_m_csf_signal_virus_EC50 = (*cfg)["nk_m_csf_signal_virus_EC50"].as<double>();
    params.nk_m_csf_signal_virus_hill = (*cfg)["nk_m_csf_signal_virus_hill"].as<double>();
    params.nk_m_csf_downreg_il_10_EC50 = (*cfg)["nk_m_csf_downreg_il_10_EC50"].as<double>();
    params.nk_m_csf_downreg_il_10_hill = (*cfg)["nk_m_csf_downreg_il_10_hill"].as<double>();

    /* infected epithelial chemokine */
    params.nk_resting_uptake_epi_inf_ck = (*cfg)["nk_resting_uptake_epi_inf_ck"].as<double>();
    params.nk_active_uptake_epi_inf_ck = (*cfg)["nk_active_uptake_epi_inf_ck"].as<double>();

    /* apoptosis chemokine */
    params.nk_secreted_apop_ck = (*cfg)["nk_secreted_apop_ck"].as<double>();

    params.nk_resting_uptake_apop_ck = (*cfg)["nk_resting_uptake_apop_ck"].as<double>();
    params.nk_active_uptake_apop_ck = (*cfg)["nk_active_uptake_apop_ck"].as<double>();

    /* phagocytosis chemokine */
    params.nk_apoptotic_secreted_phago_ck = (*cfg)["nk_apoptotic_secreted_phago_ck"].as<double>();

    /* infectious virus - no global memory; not used */

    /* non-infectious virus - no global memory; not used */

    /* movement */
    params.nk_resting_movement_rate = (*cfg)["nk_resting_movement_rate"].as<double>();
    params.nk_resting_movement_bias = (*cfg)["nk_resting_movement_bias"].as<double>();

    params.nk_active_movement_rate = (*cfg)["nk_active_movement_rate"].as<double>();
    params.nk_active_movement_bias = (*cfg)["nk_active_movement_bias"].as<double>();

    params.nk_movement_upreg_virus_EC50 = (*cfg)["nk_movement_upreg_virus_EC50"].as<double>();
    params.nk_movement_upreg_virus_hill = (*cfg)["nk_movement_upreg_virus_hill"].as<double>();

    params.nk_movement_weight_epi_inf_ck = (*cfg)["nk_movement_weight_epi_inf_ck"].as<double>();
    params.nk_movement_weight_apop_ck = (*cfg)["nk_movement_weight_apop_ck"].as<double>();

    /* recruitment */
    params.nk_recruitment_weight_il_6 = (*cfg)["nk_recruitment_weight_il_6"].as<double>();
    params.nk_recruitment_weight_il_10 = (*cfg)["nk_recruitment_weight_il_10"].as<double>();
    params.nk_recruitment_weight_il_12 = (*cfg)["nk_recruitment_weight_il_12"].as<double>();
    params.nk_recruitment_weight_epi_inf_ck = (*cfg)["nk_recruitment_weight_epi_inf_ck"].as<double>();
    params.nk_recruitment_weight_apop_ck = (*cfg)["nk_recruitment_weight_apop_ck"].as<double>();

    params.nk_recruitment_EC50 = (*cfg)["nk_recruitment_EC50"].as<double>();
    params.nk_recruitment_hill = (*cfg)["nk_recruitment_hill"].as<double>();

    /* activation and deactivation */
    params.nk_activate_weight_ifn_1 = (*cfg)["nk_activate_weight_ifn_1"].as<double>();
    params.nk_activate_weight_il_6 = (*cfg)["nk_activate_weight_il_6"].as<double>();
    params.nk_activate_weight_il_10 = (*cfg)["nk_activate_weight_il_10"].as<double>();
    params.nk_activate_weight_il_12 = (*cfg)["nk_activate_weight_il_12"].as<double>();
    params.nk_activate_weight_epi_inf_ck = (*cfg)["nk_activate_weight_epi_inf_ck"].as<double>();
    params.nk_activate_weight_apop_ck = (*cfg)["nk_activate_weight_apop_ck"].as<double>();
    params.nk_activate_weight_virus_inf = (*cfg)["nk_activate_weight_virus_inf"].as<double>();

    params.nk_activate_EC50 = (*cfg)["nk_activate_EC50"].as<double>();
    params.nk_activate_hill = (*cfg)["nk_activate_hill"].as<double>();

    /* apoptosis */
    params.nk_apoptosis_timescale = (*cfg)["nk_apoptosis_timescale"].as<double>();
    params.nk_apoptosis_max_rate = (*cfg)["nk_apoptosis_max_rate"].as<double>();
    params.nk_apoptosis_half_max = (*cfg)["nk_apoptosis_half_max"].as<double>();
    params.nk_apoptosis_hill = (*cfg)["nk_apoptosis_hill"].as<double>();
    params.nk_apoptosis_by_ifn_max_rate = (*cfg)["nk_apoptosis_by_ifn_max_rate"].as<double>();
    params.nk_apoptosis_by_ifn_half_max = (*cfg)["nk_apoptosis_by_ifn_half_max"].as<double>();
    params.nk_apoptosis_by_ifn_hill = (*cfg)["nk_apoptosis_by_ifn_hill"].as<double>();
    params.nk_active_apoptose_epithelial_probability = (*cfg)["nk_active_apoptose_epithelial_probability"].as<double>();
    params.nk_active_apoptose_immune_probability = (*cfg)["nk_active_apoptose_immune_probability"].as<double>();

    /* phagocytosis */
    params.nk_phagocytosis_timescale = (*cfg)["nk_phagocytosis_timescale"].as<double>();

  } catch (std::exception &e) {
    CAM_ERROR(e.what(), CAM_ERROR_IO);
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * NK constructor
 * =============================================================================
 */
NK::NK (const int x, const int y, const NKType tp, const NKState st,
        struct cam_rd_sys* rds, const bool use_output, int* status)
  : Agent{ x, y, AgentType::nk, rds }
{
  /* agent is constructed first */

  /* set nk type */
  set_type(tp);

  NKState state;
  if (st == NKState::null) {
    /* calculate random number between 0 and 1 */
    double r = distribution_0_1(rng);
    if (r < params.nk_active_probability) {
      /* nk state is active */
      state = NKState::active;
    } else {
      /* nk state is resting */
      state = NKState::resting;
    }
  } else {
    state = st;
  }

  /* set nk state, age, lifespan, secretions and uptake */
  set_state(state);

  /* set the area */
  double nk_area = CAM_PI * params.nk_radius * params.nk_radius;
  set_area(nk_area);

  /* update the reaction diffusion secretions & uptakes */
  inc_rds_values();

  /* resident on cell */
  resident_on_cell = nullptr;

  /* apoptosis */
  is_being_apoptosed = false;
  apoptosis_start_time = 0.0;
  apoptosis_timescale = params.nk_apoptosis_timescale;
  num_cells_apoptosing = 0;
  is_apoptosing_epithelial = false;
  apoptosing_agent = nullptr;

  /* phagocytosis */
  phagocytosis_timescale = params.nk_phagocytosis_timescale;
  is_phagocytosing_epithelial = false;
  phagocytosing_agent = nullptr;
  phagocytosed_by = nullptr;

  std::string outf;
  if (use_output) {
    outf = output_folder + std::string("/path_nk.txt");
    output = fopen(outf.c_str(), "w");

    outf = output_folder + std::string("/nk_total_virus_uptaken.txt");
    uptake_output = fopen(outf.c_str(), "w");

    outf = output_folder + std::string("/nk_lifespan.txt");
    lifespan_output = fopen(outf.c_str(), "w");

    outf = output_folder + std::string("/nk_movement_rate.txt");
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
 * NK destructor
 * =============================================================================
 */
NK::~NK ()
{
  /* need to decrease the cytokine secretions and uptakes */
  dec_rds_values();

  if (resident_on_cell != nullptr) {
    resident_on_cell = nullptr;
  }

  if (apoptosing_agent != nullptr) {
    apoptosing_agent = nullptr;
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
 * set nk state, reset age, lifespan, secretions and uptakes
 * =============================================================================
 */
void
NK::set_state (const NKState st)
{
  /* set state */
  nk_state = st;

  /* reset age */
  set_age(0.0);

  /* set the lifespan, depending on type and state */
  double lspan = 0.0;
  switch (nk_state) {
    case NKState::resting:
      lspan = params.nk_resting_lifespan;
      break;
    case NKState::active:
      lspan = params.nk_active_lifespan;
      break;
    default:
      break;
  }
  set_lifespan(lspan, true);

  /* set secretion depending on type and state */
  set_secretion(st);

  /* set uptake depending on type and state */
  set_uptake(st);

  /* set movement rates and bias */
  set_movement(st);
}


/*
 * =============================================================================
 * set the cytokine and chemokine secretion based on state
 * =============================================================================
 */
void
NK::set_secretion (const NKState st)
{
  switch (st) {
    case NKState::resting:
      basal_secreted_ifn_1 = params.nk_resting_basal_secreted_ifn_1;
      secreted_ifn_1 = params.nk_resting_secreted_ifn_1;

      basal_secreted_il_6 = params.nk_resting_basal_secreted_il_6;
      secreted_il_6 = params.nk_resting_secreted_il_6;

      basal_secreted_il_10 = params.nk_resting_basal_secreted_il_10;
      secreted_il_10 = params.nk_resting_secreted_il_10;

      basal_secreted_il_12 = params.nk_resting_basal_secreted_il_12;
      secreted_il_12 = params.nk_resting_secreted_il_12;

      basal_secreted_g_csf = params.nk_resting_basal_secreted_g_csf;
      secreted_g_csf = params.nk_resting_secreted_g_csf;

      basal_secreted_m_csf = params.nk_resting_basal_secreted_m_csf;
      secreted_m_csf = params.nk_resting_secreted_m_csf;

      basal_secreted_epi_inf_ck = 0.0;
      secreted_epi_inf_ck = 0.0;

      secreted_apop_ck = params.nk_secreted_apop_ck;
      secreted_phago_ck = 0.0;

      break;
    case NKState::active:
      basal_secreted_ifn_1 = params.nk_active_basal_secreted_ifn_1;
      secreted_ifn_1 = params.nk_active_secreted_ifn_1;

      basal_secreted_il_6 = params.nk_active_basal_secreted_il_6;
      secreted_il_6 = params.nk_active_secreted_il_6;

      basal_secreted_il_10 = params.nk_active_basal_secreted_il_10;
      secreted_il_10 = params.nk_active_secreted_il_10;

      basal_secreted_il_12 = params.nk_active_basal_secreted_il_12;
      secreted_il_12 = params.nk_active_secreted_il_12;

      basal_secreted_g_csf = params.nk_active_basal_secreted_g_csf;
      secreted_g_csf = params.nk_active_secreted_g_csf;

      basal_secreted_m_csf = params.nk_active_basal_secreted_m_csf;
      secreted_m_csf = params.nk_active_secreted_m_csf;

      basal_secreted_epi_inf_ck = 0.0;
      secreted_epi_inf_ck = 0.0;

      secreted_apop_ck = params.nk_secreted_apop_ck;
      secreted_phago_ck = 0.0;

      break;
    case NKState::apoptotic:
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

      secreted_phago_ck = params.nk_apoptotic_secreted_phago_ck;

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
NK::set_uptake (const NKState st)
{
  switch (st) {
    case NKState::resting:
      uptake_ifn_1 = params.nk_resting_uptake_ifn_1;

      uptake_il_6 = params.nk_resting_uptake_il_6;
      uptake_il_10 = params.nk_resting_uptake_il_10;
      uptake_il_12 = params.nk_resting_uptake_il_12;

      uptake_g_csf = params.nk_resting_uptake_g_csf;
      uptake_m_csf = params.nk_resting_uptake_m_csf;

      uptake_epi_inf_ck = params.nk_resting_uptake_epi_inf_ck;
      uptake_apop_ck = params.nk_resting_uptake_apop_ck;
      uptake_phago_ck = 0.0;

      break;
    case NKState::active:
      uptake_ifn_1 = params.nk_active_uptake_ifn_1;

      uptake_il_6 = params.nk_active_uptake_il_6;
      uptake_il_10 = params.nk_active_uptake_il_10;
      uptake_il_12 = params.nk_active_uptake_il_12;

      uptake_g_csf = params.nk_active_uptake_g_csf;
      uptake_m_csf = params.nk_active_uptake_m_csf;

      uptake_epi_inf_ck = params.nk_active_uptake_epi_inf_ck;
      uptake_apop_ck = params.nk_active_uptake_apop_ck;
      uptake_phago_ck = 0.0;

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

      break;
  }

  regulate_uptake();
}


void
NK::set_movement (const NKState st)
{
  switch (st) {
    case NKState::resting:
      movement_rate = params.nk_resting_movement_rate;
      movement_bias = params.nk_resting_movement_bias;
      break;
    case NKState::active:
      movement_rate = params.nk_active_movement_rate;
      movement_bias = params.nk_active_movement_bias;
      break;
    default:
      movement_rate = 0.0;
      movement_bias = 0.0;
      break;
  }
}


void
NK::change_state (const NKState st)
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
NK::update_lifespan (const NKState st)
{
  (void) st;

  return;

  // double inhib = cam_util_hill_function(
  //   params.nk_lifespan_downreg_virus_EC50, total_virus_uptaken, params.nk_lifespan_downreg_virus_hill
  // );
  //
  // double lspan = get_lifespan(true)*inhib; // get_lifespan(true) returns max lifespan
  // set_lifespan(lspan, false);
}


void
NK::update_movement (const NKState st)
{
  /* reset the movement rate */
  set_movement(st);

  return;

  // double conc = cam_util_hill_function(
  //   total_virus_uptaken, params.nk_movement_upreg_virus_EC50, params.nk_movement_upreg_virus_hill
  // );
  //
  // movement_rate *= (1.0 + conc); // NEEDSWORK: => max at 2*movement_rate
}


void
NK::regulate_secretion ()
{
  double conc = 0.0;
  double signal = 0.0;
  double upreg = 0.0;
  double downreg = 0.0;

  std::pair<int, int> xy = get_coordinates();

  const double virus_inf = rd_sys->virus_inf->get_conc(xy.first, xy.second);
  const double ifn_1 = rd_sys->ifn_1->get_conc(xy.first, xy.second);
  const double il_6 = rd_sys->il_6->get_conc(xy.first, xy.second);
  const double il_10 = rd_sys->il_10->get_conc(xy.first, xy.second);
  const double il_12 = rd_sys->il_12->get_conc(xy.first, xy.second);
  const double g_csf = rd_sys->g_csf->get_conc(xy.first, xy.second);
  const double m_csf = rd_sys->m_csf->get_conc(xy.first, xy.second);

  /* T1 IFN */
  conc = virus_inf;
  signal = cam_util_hill_function(
    conc, params.nk_ifn_1_signal_virus_EC50, params.nk_ifn_1_signal_virus_hill
  );

  conc = il_10;
  downreg = cam_util_hill_function(
    params.nk_ifn_1_downreg_il_10_EC50, conc, params.nk_ifn_1_downreg_il_10_hill
  );

  conc = ifn_1;
  upreg = cam_util_hill_function(
    conc, params.nk_ifn_1_EC50, params.nk_ifn_1_hill
  );

  secreted_ifn_1 = (
    ((basal_secreted_ifn_1 / (1.0 + conc*conc)) + secreted_ifn_1*upreg*downreg)*signal
  );


  /* IL-6 */
  conc = virus_inf;
  signal = cam_util_hill_function(
    conc, params.nk_il_6_signal_virus_EC50, params.nk_il_6_signal_virus_hill
  );

  conc = il_10;
  downreg = cam_util_hill_function(
    params.nk_il_6_downreg_il_10_EC50, conc, params.nk_il_6_downreg_il_10_hill
  );

  conc = il_6;
  upreg = cam_util_hill_function(
    conc, params.nk_il_6_EC50, params.nk_il_6_hill
  );

  secreted_il_6 = (
    ((basal_secreted_il_6 / (1.0 + conc*conc)) + secreted_il_6*upreg*downreg)*signal
  );


  /* IL-12 */
  conc = virus_inf;
  signal = cam_util_hill_function(
    conc, params.nk_il_12_signal_virus_EC50, params.nk_il_12_signal_virus_hill
  );

  conc = il_10;
  downreg = cam_util_hill_function(
    params.nk_il_12_downreg_il_10_EC50, conc, params.nk_il_12_downreg_il_10_hill
  );

  conc = il_12;
  upreg = cam_util_hill_function(
    conc, params.nk_il_12_EC50, params.nk_il_12_hill
  );

  secreted_il_12 = (
    ((basal_secreted_il_12 / (1.0 + conc*conc)) + secreted_il_12*upreg*downreg)*signal
  );


  /* IL-10 */
  conc = il_6;
  signal = cam_util_hill_function(
    conc, params.nk_il_10_signal_il_6_EC50, params.nk_il_10_signal_il_6_hill
  );
  conc = il_12;
  signal += cam_util_hill_function(
    conc, params.nk_il_10_signal_il_12_EC50, params.nk_il_10_signal_il_12_hill
  );
  conc = g_csf;
  signal += cam_util_hill_function(
    conc, params.nk_il_10_signal_g_csf_EC50, params.nk_il_10_signal_g_csf_hill
  );
  conc = m_csf;
  signal += cam_util_hill_function(
    conc, params.nk_il_10_signal_m_csf_EC50, params.nk_il_10_signal_m_csf_hill
  );

  conc = il_10;
  upreg = cam_util_hill_function(
    conc, params.nk_il_10_EC50, params.nk_il_10_hill
  );

  secreted_il_10 = (
    ((basal_secreted_il_10 / (1.0 + conc*conc)) + secreted_il_10*upreg)*signal
  );


  /* G-CSF */
  conc = virus_inf;
  signal = cam_util_hill_function(
    conc, params.nk_g_csf_signal_virus_EC50, params.nk_g_csf_signal_virus_hill
  );

  conc = il_10;
  downreg = cam_util_hill_function(
    params.nk_g_csf_downreg_il_10_EC50, conc, params.nk_g_csf_downreg_il_10_hill
  );

  conc = g_csf;
  upreg = cam_util_hill_function(
    conc, params.nk_g_csf_EC50, params.nk_g_csf_hill
  );

  secreted_g_csf = (
    ((basal_secreted_g_csf / (1.0 + conc*conc)) + secreted_g_csf*upreg*downreg)*signal
  );


  /* M-CSF */
  conc = virus_inf;
  signal = cam_util_hill_function(
    conc, params.nk_m_csf_signal_virus_EC50, params.nk_m_csf_signal_virus_hill
  );

  conc = il_10;
  downreg = cam_util_hill_function(
    params.nk_m_csf_downreg_il_10_EC50, conc, params.nk_m_csf_downreg_il_10_hill
  );

  conc = m_csf;
  upreg = cam_util_hill_function(
    conc, params.nk_m_csf_EC50, params.nk_m_csf_hill
  );

  secreted_m_csf = (
    ((basal_secreted_m_csf / (1.0 + conc*conc)) + secreted_m_csf*upreg*downreg)*signal
  );


  /* infected epithelial chemokine - not used */
}


void
NK::regulate_uptake ()
{
  /* nothing to do */
}


void
NK::update_secretions ()
{
  /* update the reaction diffusion secretions and uptakes */
  dec_rds_values();

  /* update the lifespan */
  update_lifespan(nk_state);

  /* reset the secretions */
  set_secretion(nk_state);

  /* reset the uptakes */
  set_uptake(nk_state);

  /* update the movement */
  update_movement(nk_state);

  /* update the reaction diffusion secretions and uptakes */
  inc_rds_values();
}


/*
 * =============================================================================
 * Compute the nk action (activation, deactivation, triggering apoptosis)
 * =============================================================================
 */
void
NK::action (const double time)
{
  /* random number */
  double r;

  /* coords */
  std::pair<int, int> xy = get_coordinates();

  if (nk_state == NKState::apoptotic || is_being_apoptosed || is_being_phagocytosed || nk_state == NKState::phagocytosed) {
    return;
  }

  /* update the secretions */
  update_secretions();

  double ifn_uptake = rd_sys->ifn_1->get_conc(xy.first, xy.second);
  double trigger_apop = params.nk_apoptosis_by_ifn_max_rate * cam_util_hill_function(ifn_uptake, params.nk_apoptosis_by_ifn_half_max, params.nk_apoptosis_by_ifn_hill);

  r = distribution_0_1(rng);
  if (r < trigger_apop) {
    increase_age(10.0*get_lifespan(true));

    return;
  }

  /* get the aggregate of the cytokines */
  double cyt_tmp = (
    params.nk_activate_weight_virus_inf * rd_sys->virus_inf->get_conc(xy.first, xy.second) +
    params.nk_activate_weight_ifn_1 * rd_sys->ifn_1->get_conc(xy.first, xy.second) +
    params.nk_activate_weight_il_6 * rd_sys->il_6->get_conc(xy.first, xy.second) +
    params.nk_activate_weight_il_12 * rd_sys->il_12->get_conc(xy.first, xy.second) +
    params.nk_activate_weight_epi_inf_ck * rd_sys->epi_inf_ck->get_conc(xy.first, xy.second) +
    params.nk_activate_weight_apop_ck * rd_sys->apop_ck->get_conc(xy.first, xy.second) -
    params.nk_activate_weight_il_10 * rd_sys->il_10->get_conc(xy.first, xy.second)
  );

  double cyt = cam_util_hill_function(
    cyt_tmp, params.nk_activate_EC50, params.nk_activate_hill
  );

  /* resting actions */
  if (nk_state == NKState::resting) {
    /* activation */
    r = distribution_0_1(rng);
    if (r < cyt && !is_apoptosing_epithelial && apoptosing_agent == nullptr && !is_being_apoptosed && !is_being_phagocytosed) {
      change_state(NKState::active);

      /* update the secretions */
      update_secretions();

      return;
    }
  }


  /* active actions */
  if (nk_state == NKState::active) {
    /* apoptose the epithelial cell the NK cell is resident on */
    int ret_code = apoptose_epithelial(time);
    if (ret_code == 0) {
      /* successful apoptosis action; start, continue or end */
      return;
    }

    /* NK cell has not triggered, continued or ended apoptosis in epithelial
     * cell. Therefore, trigger apoptosis in resident inflammatory immune cells */
    apoptose_immune(time);

    /* deactivation */
    r = distribution_0_1(rng);
    if (r < cyt && !is_apoptosing_epithelial && apoptosing_agent == nullptr && !is_being_apoptosed && !is_being_phagocytosed) {
      change_state(NKState::resting);

      /* update the secretions */
      update_secretions();

      return;
    }
  }
}


/*
 * =============================================================================
 * check & trigger apoptosis of epithelial cell
 *
 * return codes == 0 - successfull apoptosis action
 *              != 0 - available for immune cell apoptosis
 * =============================================================================
 */
int
NK::apoptose_epithelial (const double time)
{
  if (apoptosing_agent != nullptr) {
    /* this NK cell is already apoptosing an immune cell */
    return -1;
  }


  /* apoptosis of epithelial cell */
  if (is_apoptosing_epithelial && !resident_on_cell->is_being_apoptosed) {
    /* apoptosis is complete */

    is_apoptosing_epithelial = false;

    return 0;

  } else if (!is_apoptosing_epithelial && resident_on_cell->is_being_apoptosed) {
    /* apoptosis already started elsewhere */

    /* change time scale of apoptosis due to NK action */
    double r = distribution_0_1(rng);
    if (r < (params.nk_active_apoptose_epithelial_probability)) {
      resident_on_cell->num_cells_apoptosing++;

      double tscale = apoptosis_timescale;

      resident_on_cell->set_apoptosis_timescale(tscale);

      is_apoptosing_epithelial = true;

      return 0;
    }

  } else if (!is_apoptosing_epithelial && !resident_on_cell->is_being_apoptosed) {
    /* NK cell is starting the apoptosis */

    if (resident_on_cell->check_state(EpithelialState::infected)) {
      double r = distribution_0_1(rng);
      if (r < (params.nk_active_apoptose_epithelial_probability)) {
        /* start apoptosis */

        double tscale = apoptosis_timescale;

        resident_on_cell->start_apoptosis(time);
        resident_on_cell->set_apoptosis_timescale(tscale);
        resident_on_cell->switch_flag_apoptosis();
        is_apoptosing_epithelial = true;

        return 0;
      }
    }
  }

  /* want to check whether immune cell apoptosis can be triggered */
  return -2;
}


void
NK::apoptose_immune (const double time)
{
  if (is_apoptosing_epithelial) {
    /* this NK cell is apoptosing an epithelial cell then return */
    return;
  }


  /* apoptosis of resident immune cell */
  if (apoptosing_agent == nullptr) {
    /* starting apoptosis of immune cell */

    std::vector<Agent*>* agents = resident_on_cell->get_residents();
    for (int i = 0; i < (int) (*agents).size(); i++) {
      Agent* ag = (*agents)[i];
      if (ag == this) {
        /* NOTE: might have to dynamic_cast "this" object to Agent* */
        /* if the agent is the current object */
        continue;
      }

      if (!ag->is_active()) {
        /* only active immune cells can be apoptosed */
        continue;
      }

      if (ag->is_being_apoptosed) {
        /* apoptosis already started elsewhere */

        double r = distribution_0_1(rng);
        if (r < (params.nk_active_apoptose_immune_probability)) {
          ag->num_cells_apoptosing++;

          double tscale = apoptosis_timescale;

          /* change time scale of apoptosis due to NK action */
          ag->set_apoptosis_timescale(tscale);

          apoptosing_agent = ag;
        }

        return;
      }

      double r = distribution_0_1(rng);
      if (r < (params.nk_active_apoptose_immune_probability)) {
        /* start apoptosis */

        double tscale = apoptosis_timescale;

        ag->start_apoptosis(time);
        ag->set_apoptosis_timescale(tscale);

        apoptosing_agent = ag;

        return;
      }
    }

  } else {
    if (!apoptosing_agent->is_being_apoptosed) {
      /* end of apoptosis */

      apoptosing_agent = nullptr;
    } /* else { apoptosis ongoing } */
  }
}


int
NK::migrate (const double dt, const double tm, std::function<int(Agent*, std::vector<std::pair<int,int>>)> fndsp, std::function<bool(Agent*, std::pair<int,int>)> chksp)
{
  Epithelial* epi = resident_on_cell;
  const int sz = (int) epi->immediate_moore.size();

  /* return index */
  int ind = sz;

  /* movement rate */
  if (is_apoptosing_epithelial || apoptosing_agent != nullptr || is_being_phagocytosed ||
      is_being_apoptosed || (nk_state == NKState::apoptotic)) {
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
        params.nk_movement_weight_epi_inf_ck * rd_sys->epi_inf_ck->get_conc(xy.first, xy.second) +
        params.nk_movement_weight_apop_ck * rd_sys->apop_ck->get_conc(xy.first, xy.second)
      );

      total_cyt += pow(nbr_cyt, movement_bias);
    }

    if ((total_cyt < CAM_TOL) || nk_state == NKState::resting) {
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
          params.nk_movement_weight_epi_inf_ck * rd_sys->epi_inf_ck->get_conc(xy.first, xy.second) +
          params.nk_movement_weight_apop_ck * rd_sys->apop_ck->get_conc(xy.first, xy.second)
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
NK::recruit (std::vector<std::pair<int,int>> nbrs)
{
  if (params.nk_recruitment_hill < CAM_TOL) {
    return false;
  }

  /* Get the cytokine levels from all neighbours */
  double total_cyt = 0.0;
  double nbr_il_12 = 0.0, nbr_il_6 = 0.0, nbr_epi_inf_ck = 0.0, nbr_apop_ck = 0.0, nbr_il_10 = 0.0;

  for (std::pair<int,int> nb:nbrs) {
    nbr_il_12 += rd_sys->il_12->get_conc(nb.first, nb.second);
    nbr_il_6 += rd_sys->il_6->get_conc(nb.first, nb.second);
    nbr_epi_inf_ck += rd_sys->epi_inf_ck->get_conc(nb.first, nb.second);
    nbr_apop_ck += rd_sys->apop_ck->get_conc(nb.first, nb.second);
    nbr_il_10 += rd_sys->il_10->get_conc(nb.first, nb.second);
  }

  total_cyt = (
    params.nk_recruitment_weight_il_12 * nbr_il_12 +
    params.nk_recruitment_weight_il_6 * nbr_il_6 +
    params.nk_recruitment_weight_epi_inf_ck * nbr_epi_inf_ck +
    params.nk_recruitment_weight_apop_ck * nbr_apop_ck -
    params.nk_recruitment_weight_il_10 * nbr_il_10
  ) / ((double) nbrs.size());

  double cyt = cam_util_hill_function(
    total_cyt, params.nk_recruitment_EC50, params.nk_recruitment_hill
  );

  double r = distribution_0_1(rng);
  if (r < cyt) {
    return true;
  }

  return false;
}


/*
 * =============================================================================
 * Start the apoptosis in the NK cell
 * =============================================================================
 */
void
NK::start_apoptosis (const double time)
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
 * End the apoptosis in the NK cell
 * =============================================================================
 */
void
NK::end_apoptosis ()
{
  std::pair<int, int> xy = get_coordinates();

  is_being_apoptosed = false;
  apoptosis_start_time = 0.0;
  num_cells_apoptosing = 0;

  change_state(NKState::apoptotic);

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
NK::set_apoptosis_timescale (const double tms)
{
  double timescale = tms / ((double) num_cells_apoptosing);

  /* change time scale of apoptosis due to NK action */
  if (timescale < apoptosis_timescale) {
    apoptosis_timescale = timescale;
  }
}


/*
 * =============================================================================
 * Start the phagocytosis in the NK cell
 * =============================================================================
 */
void
NK::start_phagocytosis (const double time)
{
  is_being_phagocytosed = true;
  phagocytosis_start_time = time;
}


/*
 * =============================================================================
 * End the phagocytosis in the NK cell
 * =============================================================================
 */
void
NK::end_phagocytosis ()
{
  is_being_phagocytosed = false;
  phagocytosis_start_time = 0.0;

  phagocytosed_by->phagocytosing_agent = nullptr;
  phagocytosed_by = nullptr;
  phagocytosing_agent = nullptr;

  change_state(NKState::phagocytosed);

  /* update the secretions */
  update_secretions();
}


/*
 * =============================================================================
 * Set the phagocytosis timescale
 * =============================================================================
 */
void
NK::set_phagocytosis_timescale (const double tms)
{
  /* change time scale of phagocytosis */
  if (tms < phagocytosis_timescale) {
    phagocytosis_timescale = tms;
  }
}
