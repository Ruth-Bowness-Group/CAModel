#include <cstring>

#include "cam_epithelial.h"
#include "epithelial.h"

#include "../cam_error.h"
#include "../cam_parameters.h"
#include "../cam_random.h"
#include "../cam_util.h"


/* struct containing epithelial specific parameters */
static struct epithelial_params params;


/*
 * =============================================================================
 * Set the epithelial specific parameters
 * =============================================================================
 */
int
cam_epithelial_set_params (YAML::Node* cfg)
{
  int status = 0;
  std::string vrn, ren;

  try {
    /* spatial step (corresponds to epithelial width and height) */
    params.epithelial_width = (*cfg)["rd_spatial_step"].as<double>();
    params.epithelial_height = (*cfg)["rd_spatial_step"].as<double>();

    /*
     * cytokines
     */

    /* T1 IFN */
    params.epithelial_healthy_basal_secreted_ifn_1 = (*cfg)["epithelial_healthy_basal_secreted_ifn_1"].as<double>();
    params.epithelial_healthy_secreted_ifn_1 = (*cfg)["epithelial_healthy_secreted_ifn_1"].as<double>();
    params.epithelial_healthy_uptake_ifn_1 = (*cfg)["epithelial_healthy_uptake_ifn_1"].as<double>();

    params.epithelial_infected_basal_secreted_ifn_1 = (*cfg)["epithelial_infected_basal_secreted_ifn_1"].as<double>();
    params.epithelial_infected_secreted_ifn_1 = (*cfg)["epithelial_infected_secreted_ifn_1"].as<double>();
    params.epithelial_infected_uptake_ifn_1 = (*cfg)["epithelial_infected_uptake_ifn_1"].as<double>();

    params.epithelial_ifn_1_EC50 = (*cfg)["epithelial_ifn_1_EC50"].as<double>();
    params.epithelial_ifn_1_hill = (*cfg)["epithelial_ifn_1_hill"].as<double>();
    params.epithelial_ifn_1_signal_virus_EC50 = (*cfg)["epithelial_ifn_1_signal_virus_EC50"].as<double>();
    params.epithelial_ifn_1_signal_virus_hill = (*cfg)["epithelial_ifn_1_signal_virus_hill"].as<double>();
    params.epithelial_ifn_1_downreg_il_10_EC50 = (*cfg)["epithelial_ifn_1_downreg_il_10_EC50"].as<double>();
    params.epithelial_ifn_1_downreg_il_10_hill = (*cfg)["epithelial_ifn_1_downreg_il_10_hill"].as<double>();

    params.epithelial_infected_delay_ifn_1 = (*cfg)["epithelial_infected_delay_ifn_1"].as<double>();

    /* IL-6 */
    params.epithelial_healthy_basal_secreted_il_6 = (*cfg)["epithelial_healthy_basal_secreted_il_6"].as<double>();
    params.epithelial_healthy_secreted_il_6 = (*cfg)["epithelial_healthy_secreted_il_6"].as<double>();
    params.epithelial_healthy_uptake_il_6 = (*cfg)["epithelial_healthy_uptake_il_6"].as<double>();

    params.epithelial_infected_basal_secreted_il_6 = (*cfg)["epithelial_infected_basal_secreted_il_6"].as<double>();
    params.epithelial_infected_secreted_il_6 = (*cfg)["epithelial_infected_secreted_il_6"].as<double>();
    params.epithelial_infected_uptake_il_6 = (*cfg)["epithelial_infected_uptake_il_6"].as<double>();

    params.epithelial_il_6_EC50 = (*cfg)["epithelial_il_6_EC50"].as<double>();
    params.epithelial_il_6_hill = (*cfg)["epithelial_il_6_hill"].as<double>();
    params.epithelial_il_6_signal_virus_EC50 = (*cfg)["epithelial_il_6_signal_virus_EC50"].as<double>();
    params.epithelial_il_6_signal_virus_hill = (*cfg)["epithelial_il_6_signal_virus_hill"].as<double>();
    params.epithelial_il_6_downreg_il_10_EC50 = (*cfg)["epithelial_il_6_downreg_il_10_EC50"].as<double>();
    params.epithelial_il_6_downreg_il_10_hill = (*cfg)["epithelial_il_6_downreg_il_10_hill"].as<double>();

    /* IL-10 */
    params.epithelial_healthy_basal_secreted_il_10 = (*cfg)["epithelial_healthy_basal_secreted_il_10"].as<double>();
    params.epithelial_healthy_secreted_il_10 = (*cfg)["epithelial_healthy_secreted_il_10"].as<double>();
    params.epithelial_healthy_uptake_il_10 = (*cfg)["epithelial_healthy_uptake_il_10"].as<double>();

    params.epithelial_infected_basal_secreted_il_10 = (*cfg)["epithelial_infected_basal_secreted_il_10"].as<double>();
    params.epithelial_infected_secreted_il_10 = (*cfg)["epithelial_infected_secreted_il_10"].as<double>();
    params.epithelial_infected_uptake_il_10 = (*cfg)["epithelial_infected_uptake_il_10"].as<double>();

    params.epithelial_il_10_EC50 = (*cfg)["epithelial_il_10_EC50"].as<double>();
    params.epithelial_il_10_hill = (*cfg)["epithelial_il_10_hill"].as<double>();
    params.epithelial_il_10_signal_il_6_EC50 = (*cfg)["epithelial_il_10_signal_il_6_EC50"].as<double>();
    params.epithelial_il_10_signal_il_6_hill = (*cfg)["epithelial_il_10_signal_il_6_hill"].as<double>();
    params.epithelial_il_10_signal_il_12_EC50 = (*cfg)["epithelial_il_10_signal_il_12_EC50"].as<double>();
    params.epithelial_il_10_signal_il_12_hill = (*cfg)["epithelial_il_10_signal_il_12_hill"].as<double>();
    params.epithelial_il_10_signal_g_csf_EC50 = (*cfg)["epithelial_il_10_signal_g_csf_EC50"].as<double>();
    params.epithelial_il_10_signal_g_csf_hill = (*cfg)["epithelial_il_10_signal_g_csf_hill"].as<double>();
    params.epithelial_il_10_signal_m_csf_EC50 = (*cfg)["epithelial_il_10_signal_m_csf_EC50"].as<double>();
    params.epithelial_il_10_signal_m_csf_hill = (*cfg)["epithelial_il_10_signal_m_csf_hill"].as<double>();

    /* IL-12 */
    params.epithelial_healthy_basal_secreted_il_12 = (*cfg)["epithelial_healthy_basal_secreted_il_12"].as<double>();
    params.epithelial_healthy_secreted_il_12 = (*cfg)["epithelial_healthy_secreted_il_12"].as<double>();
    params.epithelial_healthy_uptake_il_12 = (*cfg)["epithelial_healthy_uptake_il_12"].as<double>();

    params.epithelial_infected_basal_secreted_il_12 = (*cfg)["epithelial_infected_basal_secreted_il_12"].as<double>();
    params.epithelial_infected_secreted_il_12 = (*cfg)["epithelial_infected_secreted_il_12"].as<double>();
    params.epithelial_infected_uptake_il_12 = (*cfg)["epithelial_infected_uptake_il_12"].as<double>();

    params.epithelial_il_12_EC50 = (*cfg)["epithelial_il_12_EC50"].as<double>();
    params.epithelial_il_12_hill = (*cfg)["epithelial_il_12_hill"].as<double>();
    params.epithelial_il_12_signal_virus_EC50 = (*cfg)["epithelial_il_12_signal_virus_EC50"].as<double>();
    params.epithelial_il_12_signal_virus_hill = (*cfg)["epithelial_il_12_signal_virus_hill"].as<double>();
    params.epithelial_il_12_downreg_il_10_EC50 = (*cfg)["epithelial_il_12_downreg_il_10_EC50"].as<double>();
    params.epithelial_il_12_downreg_il_10_hill = (*cfg)["epithelial_il_12_downreg_il_10_hill"].as<double>();

    /* G-CSF */
    params.epithelial_healthy_basal_secreted_g_csf = (*cfg)["epithelial_healthy_basal_secreted_g_csf"].as<double>();
    params.epithelial_healthy_secreted_g_csf = (*cfg)["epithelial_healthy_secreted_g_csf"].as<double>();
    params.epithelial_healthy_uptake_g_csf = (*cfg)["epithelial_healthy_uptake_g_csf"].as<double>();

    params.epithelial_infected_basal_secreted_g_csf = (*cfg)["epithelial_infected_basal_secreted_g_csf"].as<double>();
    params.epithelial_infected_secreted_g_csf = (*cfg)["epithelial_infected_secreted_g_csf"].as<double>();
    params.epithelial_infected_uptake_g_csf = (*cfg)["epithelial_infected_uptake_g_csf"].as<double>();

    params.epithelial_g_csf_EC50 = (*cfg)["epithelial_g_csf_EC50"].as<double>();
    params.epithelial_g_csf_hill = (*cfg)["epithelial_g_csf_hill"].as<double>();
    params.epithelial_g_csf_signal_virus_EC50 = (*cfg)["epithelial_g_csf_signal_virus_EC50"].as<double>();
    params.epithelial_g_csf_signal_virus_hill = (*cfg)["epithelial_g_csf_signal_virus_hill"].as<double>();
    params.epithelial_g_csf_downreg_il_10_EC50 = (*cfg)["epithelial_g_csf_downreg_il_10_EC50"].as<double>();
    params.epithelial_g_csf_downreg_il_10_hill = (*cfg)["epithelial_g_csf_downreg_il_10_hill"].as<double>();

    /* M-CSF */
    params.epithelial_healthy_basal_secreted_m_csf = (*cfg)["epithelial_healthy_basal_secreted_m_csf"].as<double>();
    params.epithelial_healthy_secreted_m_csf = (*cfg)["epithelial_healthy_secreted_m_csf"].as<double>();
    params.epithelial_healthy_uptake_m_csf = (*cfg)["epithelial_healthy_uptake_m_csf"].as<double>();

    params.epithelial_infected_basal_secreted_m_csf = (*cfg)["epithelial_infected_basal_secreted_m_csf"].as<double>();
    params.epithelial_infected_secreted_m_csf = (*cfg)["epithelial_infected_secreted_m_csf"].as<double>();
    params.epithelial_infected_uptake_m_csf = (*cfg)["epithelial_infected_uptake_m_csf"].as<double>();

    params.epithelial_m_csf_EC50 = (*cfg)["epithelial_m_csf_EC50"].as<double>();
    params.epithelial_m_csf_hill = (*cfg)["epithelial_m_csf_hill"].as<double>();
    params.epithelial_m_csf_signal_virus_EC50 = (*cfg)["epithelial_m_csf_signal_virus_EC50"].as<double>();
    params.epithelial_m_csf_signal_virus_hill = (*cfg)["epithelial_m_csf_signal_virus_hill"].as<double>();
    params.epithelial_m_csf_downreg_il_10_EC50 = (*cfg)["epithelial_m_csf_downreg_il_10_EC50"].as<double>();
    params.epithelial_m_csf_downreg_il_10_hill = (*cfg)["epithelial_m_csf_downreg_il_10_hill"].as<double>();

    /* infected epithelial chemokine */
    params.epithelial_infected_basal_secreted_epi_inf_ck = (*cfg)["epithelial_infected_basal_secreted_epi_inf_ck"].as<double>();
    params.epithelial_infected_secreted_epi_inf_ck = (*cfg)["epithelial_infected_secreted_epi_inf_ck"].as<double>();
    params.epithelial_infected_uptake_epi_inf_ck = (*cfg)["epithelial_infected_uptake_epi_inf_ck"].as<double>();

    params.epithelial_epi_inf_ck_EC50 = (*cfg)["epithelial_epi_inf_ck_EC50"].as<double>();
    params.epithelial_epi_inf_ck_hill = (*cfg)["epithelial_epi_inf_ck_hill"].as<double>();
    params.epithelial_epi_inf_ck_signal_virus_EC50 = (*cfg)["epithelial_epi_inf_ck_signal_virus_EC50"].as<double>();
    params.epithelial_epi_inf_ck_signal_virus_hill = (*cfg)["epithelial_epi_inf_ck_signal_virus_hill"].as<double>();

    /* apoptosis chemokine */
    params.epithelial_secreted_apop_ck = (*cfg)["epithelial_secreted_apop_ck"].as<double>();

    /* phagocytosis chemokine */
    params.epithelial_apoptotic_secreted_phago_ck = (*cfg)["epithelial_apoptotic_secreted_phago_ck"].as<double>();

    /* infectious virus - no global memory; defined in replicate */

    /* non-infectious virus - no global memory; defined in replicate */

    /* thresholds */
    params.epithelial_burst_threshold_by_virions = (*cfg)["epithelial_burst_threshold_by_virions"].as<double>();

    /* apoptosis */
    params.epithelial_apoptosis_timescale = (*cfg)["epithelial_apoptosis_timescale"].as<double>();
    params.epithelial_apoptosis_max_rate = (*cfg)["epithelial_apoptosis_max_rate"].as<double>();
    params.epithelial_apoptosis_half_max = (*cfg)["epithelial_apoptosis_half_max"].as<double>();
    params.epithelial_apoptosis_hill = (*cfg)["epithelial_apoptosis_hill"].as<double>();
    params.epithelial_apoptosis_by_ifn_max_rate = (*cfg)["epithelial_apoptosis_by_ifn_max_rate"].as<double>();
    params.epithelial_apoptosis_by_ifn_half_max = (*cfg)["epithelial_apoptosis_by_ifn_half_max"].as<double>();
    params.epithelial_apoptosis_by_ifn_hill = (*cfg)["epithelial_apoptosis_by_ifn_hill"].as<double>();

    /* phagocytosis */
    params.epithelial_phagocytosis_timescale = (*cfg)["epithelial_phagocytosis_timescale"].as<double>();

    /* infectious virus proportion */
    params.epithelial_infectious_virus_proportion = (*cfg)["epithelial_infectious_virus_proportion"].as<double>();

    /* name of the viral replication model */
    vrn = (*cfg)["epithelial_viral_replication_model"].as<std::string>();

    /* name of the receptor expression model */
    ren = (*cfg)["epithelial_receptor_expression_model"].as<std::string>();

  } catch (std::exception &e) {
    CAM_ERROR(e.what(), CAM_ERROR_IO);
  }


  /* set the requested receptor expression type */
  if (!strcmp(ren.c_str(), "pcv2_discon")) {
    /* PhysiCell version 2 (https://doi.org/10.1101/2020.04.02.019075) with discontinuous entry */
    params.re_typ = cam_re_type_pcv2_discon;

  } else if (!strcmp(ren.c_str(), "null")) {
    /* No model */
    params.re_typ = nullptr;

  } else {
    CAM_ERROR("unknown receptor expression model", CAM_ERROR_VALUE);
  }


  /* set the requested viral replication type */
  if (!strcmp(vrn.c_str(), "pcv4_discon")) {
    /* PhysiCell version 4 (https://doi.org/10.1101/2020.04.02.019075) with discontinuous export */
    params.vr_typ = cam_vr_type_pcv4_discon;

  } else if (!strcmp(ren.c_str(), "null")) {
    /* No model */
    params.vr_typ = nullptr;

  } else {
    CAM_ERROR("unknown viral replication model", CAM_ERROR_VALUE);
  }


  /* set the receptor expression specific parameters */
  status = cam_recepexpr_set_params(cfg, params.re_typ);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);


  /* set the viral replication specific parameters */
  status = cam_viralrepl_set_params(cfg, params.vr_typ);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Epithelial constructor
 * =============================================================================
 */
Epithelial::Epithelial (const int x, const int y, const EpithelialType tp, const EpithelialState st, struct cam_rd_sys* rds, int* status)
{
  /* coordinates */
  set_coordinates(x, y);

  /* associate rd_sys */
  rd_sys = rds;

  /* set epithelial cell type */
  set_type(tp);

  /* area */
  area = 0.0;
  max_area = params.epithelial_width * params.epithelial_height;

  /* initialise the receptor expression model */
  if (params.re_typ == nullptr) {
    re = nullptr;
  } else {
    bool use_output = ( ((x == grid_size/2) && (y == grid_size/2)) ? true : false );

    re = new RecepExpr(params.re_typ, use_output, status);
    CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);
  }

  /* initialise the viral replication model */
  if (params.vr_typ == nullptr) {
    vr = nullptr;
  } else {
    bool use_output = ( ((x == grid_size/2) && (y == grid_size/2)) ? true : false );

    vr = new ViralRepl(params.vr_typ, use_output, status);
    CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);
  }

  /* set epithelial cell state, age and secretions */
  epithelial_infected_delay_ifn_1 = params.epithelial_infected_delay_ifn_1;

  set_state(st);

  /* update the reaction diffusion secretions and uptakes */
  inc_rds_values();

  /* initialise flags */
  flag_bursting = false;
  flag_apoptosis = false;
  flag_phagocytosis = false;

  /* apoptosis */
  is_being_apoptosed = false;
  apoptosis_start_time = 0.0;
  apoptosis_timescale = params.epithelial_apoptosis_timescale;
  num_cells_apoptosing = 0;

  /* phagocytosis */
  is_being_phagocytosed = false;
  phagocytosis_start_time = 0.0;
  phagocytosis_timescale = params.epithelial_phagocytosis_timescale;

  /* successful */
  (*status) = CAM_SUCCESS;
}


/*
 * =============================================================================
 * Epithelial destructor
 * =============================================================================
 */
Epithelial::~Epithelial ()
{
  int i = 0;

  for (i = 0; i < (int) residents.size(); i++) {
    residents[i] = nullptr;
  }

  if (re != nullptr) {
    delete re;
  }

  if (vr != nullptr) {
    delete vr;
  }

  if (rd_sys != nullptr) {
    rd_sys = nullptr;
  }
}


/*
 * =============================================================================
 * set epithelial state, reset age and set secretion of cytokines & chemokines
 * =============================================================================
 */
void
Epithelial::set_state (const EpithelialState st)
{
  /* set state */
  state = st;

  /* reset age */
  set_age(0.0);

  /* set secretion depending on state */
  set_secretion(st);

  /* set uptake depending on state */
  set_uptake(st);
}


/*
 * =============================================================================
 * Change the state of the epithelial cell
 * =============================================================================
 */
void
Epithelial::change_state (const EpithelialState st)
{
  /* update the reaction diffusion secretions and uptakes */
  dec_rds_values();

  /* set the new state (also resets the age and cytokines secretions & uptakes) */
  set_state(st);

  /* update the reaction diffusion secretions and uptakes */
  inc_rds_values();
}


/*
 * =============================================================================
 * set the cytokine and chemokine secretion based on state
 * =============================================================================
 */
void
Epithelial::set_secretion (const EpithelialState st)
{
  switch (st) {
    case EpithelialState::healthy:
      basal_secreted_ifn_1 = params.epithelial_healthy_basal_secreted_ifn_1;
      secreted_ifn_1 = params.epithelial_healthy_secreted_ifn_1;

      basal_secreted_il_6 = params.epithelial_healthy_basal_secreted_il_6;
      secreted_il_6 = params.epithelial_healthy_secreted_il_6;

      basal_secreted_il_10 = params.epithelial_healthy_basal_secreted_il_10;
      secreted_il_10 = params.epithelial_healthy_secreted_il_10;

      basal_secreted_il_12 = params.epithelial_healthy_basal_secreted_il_12;
      secreted_il_12 = params.epithelial_healthy_secreted_il_12;

      basal_secreted_g_csf = params.epithelial_healthy_basal_secreted_g_csf;
      secreted_g_csf = params.epithelial_healthy_secreted_g_csf;

      basal_secreted_m_csf = params.epithelial_healthy_basal_secreted_m_csf;
      secreted_m_csf = params.epithelial_healthy_secreted_m_csf;

      basal_secreted_epi_inf_ck = 0.0;
      secreted_epi_inf_ck = 0.0;

      secreted_apop_ck = params.epithelial_secreted_apop_ck;
      secreted_phago_ck = 0.0;

      break;

    case EpithelialState::eclipse:
      basal_secreted_ifn_1 = ( (age > epithelial_infected_delay_ifn_1) ? params.epithelial_infected_basal_secreted_ifn_1 : 0.0 );
      secreted_ifn_1 = ( (age > epithelial_infected_delay_ifn_1) ? params.epithelial_infected_secreted_ifn_1 : 0.0 );

      basal_secreted_il_6 = params.epithelial_infected_basal_secreted_il_6;
      secreted_il_6 = params.epithelial_infected_secreted_il_6;

      basal_secreted_il_10 = params.epithelial_infected_basal_secreted_il_10;
      secreted_il_10 = params.epithelial_infected_secreted_il_10;

      basal_secreted_il_12 = params.epithelial_infected_basal_secreted_il_12;
      secreted_il_12 = params.epithelial_infected_secreted_il_12;

      basal_secreted_g_csf = params.epithelial_infected_basal_secreted_g_csf;
      secreted_g_csf = params.epithelial_infected_secreted_g_csf;

      basal_secreted_m_csf = params.epithelial_infected_basal_secreted_m_csf;
      secreted_m_csf = params.epithelial_infected_secreted_m_csf;

      basal_secreted_epi_inf_ck = params.epithelial_infected_basal_secreted_epi_inf_ck;
      secreted_epi_inf_ck = params.epithelial_infected_secreted_epi_inf_ck;

      secreted_apop_ck = params.epithelial_secreted_apop_ck;
      secreted_phago_ck = 0.0;

      break;

    case EpithelialState::infected:
      basal_secreted_ifn_1 = ( (age > epithelial_infected_delay_ifn_1) ? params.epithelial_infected_basal_secreted_ifn_1 : 0.0 );
      secreted_ifn_1 = ( (age > epithelial_infected_delay_ifn_1) ? params.epithelial_infected_secreted_ifn_1 : 0.0 );

      basal_secreted_il_6 = params.epithelial_infected_basal_secreted_il_6;
      secreted_il_6 = params.epithelial_infected_secreted_il_6;

      basal_secreted_il_10 = params.epithelial_infected_basal_secreted_il_10;
      secreted_il_10 = params.epithelial_infected_secreted_il_10;

      basal_secreted_il_12 = params.epithelial_infected_basal_secreted_il_12;
      secreted_il_12 = params.epithelial_infected_secreted_il_12;

      basal_secreted_g_csf = params.epithelial_infected_basal_secreted_g_csf;
      secreted_g_csf = params.epithelial_infected_secreted_g_csf;

      basal_secreted_m_csf = params.epithelial_infected_basal_secreted_m_csf;
      secreted_m_csf = params.epithelial_infected_secreted_m_csf;

      basal_secreted_epi_inf_ck = params.epithelial_infected_basal_secreted_epi_inf_ck;
      secreted_epi_inf_ck = params.epithelial_infected_secreted_epi_inf_ck;

      secreted_apop_ck = params.epithelial_secreted_apop_ck;
      secreted_phago_ck = 0.0;

      break;

    case EpithelialState::apoptotic:
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
      secreted_phago_ck = params.epithelial_apoptotic_secreted_phago_ck;

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
Epithelial::set_uptake (const EpithelialState st)
{
  switch (st) {
    case EpithelialState::healthy:
      uptake_ifn_1 = params.epithelial_healthy_uptake_ifn_1;

      uptake_il_6 = params.epithelial_healthy_uptake_il_6;
      uptake_il_10 = params.epithelial_healthy_uptake_il_10;
      uptake_il_12 = params.epithelial_healthy_uptake_il_12;

      uptake_g_csf = params.epithelial_healthy_uptake_g_csf;
      uptake_m_csf = params.epithelial_healthy_uptake_m_csf;

      uptake_epi_inf_ck = 0.0;
      uptake_apop_ck = 0.0;
      uptake_phago_ck = 0.0;

      break;

    case EpithelialState::eclipse:
      uptake_ifn_1 = params.epithelial_infected_uptake_ifn_1;

      uptake_il_6 = params.epithelial_infected_uptake_il_6;
      uptake_il_10 = params.epithelial_infected_uptake_il_10;
      uptake_il_12 = params.epithelial_infected_uptake_il_12;

      uptake_g_csf = params.epithelial_infected_uptake_g_csf;
      uptake_m_csf = params.epithelial_infected_uptake_m_csf;

      uptake_epi_inf_ck = params.epithelial_infected_uptake_epi_inf_ck;
      uptake_apop_ck = 0.0;
      uptake_phago_ck = 0.0;

      break;

    case EpithelialState::infected:
      uptake_ifn_1 = params.epithelial_infected_uptake_ifn_1;

      uptake_il_6 = params.epithelial_infected_uptake_il_6;
      uptake_il_10 = params.epithelial_infected_uptake_il_10;
      uptake_il_12 = params.epithelial_infected_uptake_il_12;

      uptake_g_csf = params.epithelial_infected_uptake_g_csf;
      uptake_m_csf = params.epithelial_infected_uptake_m_csf;

      uptake_epi_inf_ck = params.epithelial_infected_uptake_epi_inf_ck;
      uptake_apop_ck = 0.0;
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


/*
 * =============================================================================
 * Add resident to epithelial cell
 * =============================================================================
 */
int
Epithelial::add_resident (Agent* res)
{
  if (res == nullptr) {
    CAM_ERROR("cannot add null resident", CAM_ERROR_VALUE);
  }

  /* Attach the resident */
  residents.push_back(res);

  /* Update the residents internal position */
  residents[residents.size()-1]->set_coordinates(x, y);

  /* increase area of residents */
  area += res->get_area();

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Check whether cell has residents of particular agent type
 * =============================================================================
 */
bool
Epithelial::has_residents (AgentType tp)
{
  for (int i = 0; i < (int) residents.size(); i++) {
    if (residents[i]->get_agent_type() == tp) {
      return true;
    }
  }

  return false;
}


bool
Epithelial::has_residents_inflammatory ()
{
  for (int i = 0; i < (int) residents.size(); i++) {
    if (residents[i]->is_inflammatory()) {
      return true;
    }
  }

  return false;
}


bool
Epithelial::has_residents_apoptotic ()
{
  for (int i = 0; i < (int) residents.size(); i++) {
    if (residents[i]->is_apoptotic()) {
      return true;
    }
  }

  return false;
}


/*
 * =============================================================================
 * Get residents of particular agent type
 * =============================================================================
 */
void
Epithelial::get_residents (AgentType tp, std::vector<Agent*>* res)
{
  for (int i = 0; i < (int) residents.size(); i++) {
    if (residents[i]->get_agent_type() == tp) {
      (*res).push_back(residents[i]);
    }
  }
}


/*
 * =============================================================================
 * Remove agent from residents
 * =============================================================================
 */
int
Epithelial::remove_resident (Agent* res)
{
  if (residents.empty()) {
    CAM_ERROR("No residents to remove", CAM_ERROR_VALUE);
  }

  const int sz = (int) residents.size();

  /* Remove from residents vector */
  int del_i = 0;
  while (residents[del_i] != res) {
    del_i++;
    if (del_i >= sz) {
      CAM_ERROR("Cannot remove agent as not in residents", CAM_ERROR_VALUE);
    }
  }

  /* Swap del_i entry with last entry in vector */
  residents[del_i] = residents[sz-1];

  /* Now remove the last element */
  residents.pop_back();

  /* decrease resident area */
  area -= res->get_area();

  return CAM_SUCCESS;
}


void
Epithelial::regulate_secretion ()
{
  double conc = 0.0;
  double signal = 0.0;
  double upreg = 0.0;
  double downreg = 0.0;

  const double assembled_virions = vr->get_data();
  const double il_10 = rd_sys->il_10->get_conc(x, y);
  const double ifn_1 = rd_sys->ifn_1->get_conc(x, y);
  const double il_6 = rd_sys->il_6->get_conc(x, y);
  const double il_12 = rd_sys->il_12->get_conc(x, y);
  const double g_csf = rd_sys->g_csf->get_conc(x,y);
  const double m_csf = rd_sys->m_csf->get_conc(x,y);
  const double epi_inf_ck = rd_sys->epi_inf_ck->get_conc(x, y);

  /* T1 IFN */
  signal = cam_util_hill_function(
    assembled_virions, params.epithelial_ifn_1_signal_virus_EC50, params.epithelial_ifn_1_signal_virus_hill
  );

  downreg = cam_util_hill_function(
    params.epithelial_ifn_1_downreg_il_10_EC50, il_10, params.epithelial_ifn_1_downreg_il_10_hill
  );

  conc = ifn_1;
  upreg = cam_util_hill_function(
    ifn_1, params.epithelial_ifn_1_EC50, params.epithelial_ifn_1_hill
  );

  secreted_ifn_1 = (
    ((basal_secreted_ifn_1 / (1.0 + conc*conc)) + secreted_ifn_1*upreg*downreg)*signal
  );


  /* IL-6 */
  signal = cam_util_hill_function(
    assembled_virions, params.epithelial_il_6_signal_virus_EC50, params.epithelial_il_6_signal_virus_hill
  );

  downreg = cam_util_hill_function(
    params.epithelial_il_6_downreg_il_10_EC50, il_10, params.epithelial_il_6_downreg_il_10_hill
  );

  conc = il_6;
  upreg = cam_util_hill_function(
    il_6, params.epithelial_il_6_EC50, params.epithelial_il_6_hill
  );

  secreted_il_6 = (
    ((basal_secreted_il_6 / (1.0 + conc*conc)) + secreted_il_6*upreg*downreg)*signal
  );


  /* IL-12 */
  signal = cam_util_hill_function(
    assembled_virions, params.epithelial_il_12_signal_virus_EC50, params.epithelial_il_12_signal_virus_hill
  );

  downreg = cam_util_hill_function(
    params.epithelial_il_12_downreg_il_10_EC50, il_10, params.epithelial_il_12_downreg_il_10_hill
  );

  conc = il_12;
  upreg = cam_util_hill_function(
    il_12, params.epithelial_il_12_EC50, params.epithelial_il_12_hill
  );

  secreted_il_12 = (
    ((basal_secreted_il_12 / (1.0 + conc*conc)) + secreted_il_12*upreg*downreg)*signal
  );


  /* IL-10 */
  signal = 0.25 * cam_util_hill_function(
    il_6, params.epithelial_il_10_signal_il_6_EC50, params.epithelial_il_10_signal_il_6_hill
  );
  signal += 0.25 * cam_util_hill_function(
    il_12, params.epithelial_il_10_signal_il_12_EC50, params.epithelial_il_10_signal_il_12_hill
  );
  signal += 0.25 * cam_util_hill_function(
    g_csf, params.epithelial_il_10_signal_g_csf_EC50, params.epithelial_il_10_signal_g_csf_hill
  );
  signal += 0.25 * cam_util_hill_function(
    m_csf, params.epithelial_il_10_signal_m_csf_EC50, params.epithelial_il_10_signal_m_csf_hill
  );

  conc = il_10;
  upreg = cam_util_hill_function(
    il_10, params.epithelial_il_10_EC50, params.epithelial_il_10_hill
  );

  secreted_il_10 = (
    ((basal_secreted_il_10 / (1.0 + conc*conc)) + secreted_il_10*upreg)*signal
  );


  /* G-CSF */
  signal = cam_util_hill_function(
    assembled_virions, params.epithelial_g_csf_signal_virus_EC50, params.epithelial_g_csf_signal_virus_hill
  );

  downreg = cam_util_hill_function(
    params.epithelial_g_csf_downreg_il_10_EC50, il_10, params.epithelial_g_csf_downreg_il_10_hill
  );

  conc = g_csf;
  upreg = cam_util_hill_function(
    g_csf, params.epithelial_g_csf_EC50, params.epithelial_g_csf_hill
  );

  secreted_g_csf = (
    ((basal_secreted_g_csf / (1.0 + conc*conc)) + secreted_g_csf*upreg*downreg)*signal
  );


  /* M-CSF */
  signal = cam_util_hill_function(
    assembled_virions, params.epithelial_m_csf_signal_virus_EC50, params.epithelial_m_csf_signal_virus_hill
  );

  downreg = cam_util_hill_function(
    params.epithelial_m_csf_downreg_il_10_EC50, il_10, params.epithelial_m_csf_downreg_il_10_hill
  );

  conc = m_csf;
  upreg = cam_util_hill_function(
    m_csf, params.epithelial_m_csf_EC50, params.epithelial_m_csf_hill
  );

  secreted_m_csf = (
    ((basal_secreted_m_csf / (1.0 + conc*conc)) + secreted_m_csf*upreg*downreg)*signal
  );


  /* infected epithelial chemokine */
  signal = cam_util_hill_function(
    assembled_virions, params.epithelial_epi_inf_ck_signal_virus_EC50, params.epithelial_epi_inf_ck_signal_virus_hill
  );

  conc = epi_inf_ck;
  upreg = cam_util_hill_function(
    epi_inf_ck, params.epithelial_epi_inf_ck_EC50, params.epithelial_epi_inf_ck_hill
  );

  secreted_epi_inf_ck = (
    ((basal_secreted_epi_inf_ck / (1.0 + conc*conc)) + secreted_epi_inf_ck*upreg)*signal
  );
}


void
Epithelial::regulate_uptake ()
{
  /* nothing to do */
}


void
Epithelial::update_secretions ()
{
  /* update the reaction diffusion secretions and uptakes */
  dec_rds_values();

  /* reset the secretions */
  set_secretion(state);

  /* reset the uptakes */
  set_uptake(state);

  /* update the reaction diffusion secretions and uptakes */
  inc_rds_values();
}


/*
 * =============================================================================
 * Replication
 * =============================================================================
 */
int
Epithelial::replicate (double nvirs, double time_n, double time)
{
  int status = CAM_ERROR_MISC;

  /* return if the cell has burst */
  if (state == EpithelialState::burst ||
      state == EpithelialState::phagocytosed || is_being_phagocytosed ||
      state == EpithelialState::apoptotic || is_being_apoptosed) {
    return CAM_SUCCESS;
  }

  /* update the secretions */
  update_secretions();

  double ifn_uptake = rd_sys->ifn_1->get_conc(x, y);
  double trigger_apop = params.epithelial_apoptosis_by_ifn_max_rate * cam_util_hill_function(ifn_uptake, params.epithelial_apoptosis_by_ifn_half_max, params.epithelial_apoptosis_by_ifn_hill);

  double r = distribution_0_1(rng);
  if (r < trigger_apop) {
    flag_bursting = false;
    flag_apoptosis = true;
    flag_phagocytosis = false;

    return CAM_SUCCESS;
  }

  if (nvirs < CAM_TOL && state == EpithelialState::healthy) {
    return CAM_SUCCESS;
  }

  double num_virs = nvirs;
  double uptake = 0.0, source = 0.0;

  /* update the receptor model */
  if (re != nullptr) {
    /* set the input */
    re->set_input(&num_virs, &ifn_uptake);

    /* solve the receptor expression model */
    status = re->solve_system(time_n, time, &uptake, &source);
    CAM_ERROR_CHECK(status, CAM_SUCCESS);

    /* reset the uptake */
    rd_sys->virus_inf->set_uptake(x, y, uptake, true);

    /* epithelial state is eclipse when R_ib > 0 */
    double internal_bound = ( (state == EpithelialState::healthy) ? re->get_data() : 0.0 );
    if (internal_bound > CAM_TOL) {
      change_state(EpithelialState::eclipse);

      /* update the secretions */
      update_secretions();

      flag_bursting = false;
      flag_apoptosis = false;
      flag_phagocytosis = false;
    }
  }

  /* update the viral replication */
  if (vr != nullptr) {
    /* set the input */
    vr->set_input(&source, &ifn_uptake);

    /* solve the viral replication model */
    status = vr->solve_system(time_n, time, &source);
    CAM_ERROR_CHECK(status, CAM_SUCCESS);

    /* reset the source */
    source /= max_area;
    rd_sys->virus_inf->set_source(x, y, source*params.epithelial_infectious_virus_proportion, false);
    rd_sys->virus_ninf->set_source(x, y, source*(1.0 - params.epithelial_infectious_virus_proportion), false);

    /* get internalised virus */
    double internalised = ( (state == EpithelialState::infected) ? 0.0 : vr->get_data() );
    bool exporting = ( (state == EpithelialState::infected) ? false : vr->get_param() );
    if (exporting || (internalised >= 1.0)) {
      change_state(EpithelialState::infected);

      /* update the secretions */
      update_secretions();

      flag_bursting = false;
      flag_apoptosis = false;
      flag_phagocytosis = false;
    }

    internalised = ( (state == EpithelialState::infected) ? vr->get_data() : 0.0 );
    trigger_apop = params.epithelial_apoptosis_max_rate * cam_util_hill_function(internalised, params.epithelial_apoptosis_half_max, params.epithelial_apoptosis_hill);

    r = distribution_0_1(rng);
    if (r < trigger_apop) {
      flag_bursting = false;
      flag_apoptosis = true;
      flag_phagocytosis = false;

    } else if (internalised >= params.epithelial_burst_threshold_by_virions) { /* burst < internalised */
      flag_bursting = true;
      flag_apoptosis = false;
      flag_phagocytosis = false;
    }
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Get the intracellular components
 * =============================================================================
 */
double
Epithelial::get_interior ()
{
  if (vr == nullptr || state == EpithelialState::burst ||
      state == EpithelialState::phagocytosed || is_being_phagocytosed ||
      state == EpithelialState::apoptotic || is_being_apoptosed) {
    return 0.0;
  }

  double internalised = vr->get_data();

  return internalised;
}


/*
 * =============================================================================
 * Start the apoptosis in the epithelial cell
 * =============================================================================
 */
void
Epithelial::start_apoptosis (const double time)
{
  is_being_apoptosed = true;
  apoptosis_start_time = time;
  num_cells_apoptosing = 1;

  /* secrete the chemokine which recruits NK cells to aid in apoptosis */
  rd_sys->apop_ck->inc_source(x, y, secreted_apop_ck, true);
}


/*
 * =============================================================================
 * End the apoptosis in the epithelial cell
 * =============================================================================
 */
void
Epithelial::end_apoptosis ()
{
  is_being_apoptosed = false;
  apoptosis_start_time = 0.0;
  num_cells_apoptosing = 0;

  change_state(EpithelialState::apoptotic);

  /* update the secretions */
  update_secretions();

  flag_bursting = false;
  flag_apoptosis = false;

  /* reset the NK chemokine */
  rd_sys->apop_ck->dec_source(x, y, secreted_apop_ck, true);
}


/*
 * =============================================================================
 * Set the apoptosis timescale
 * =============================================================================
 */
void
Epithelial::set_apoptosis_timescale (const double tms)
{
  double timescale = tms / ((double) num_cells_apoptosing);

  /* change time scale of apoptosis due to NK action */
  if (CAM_CMP(timescale, apoptosis_timescale, CAM_TOL)) { /* tms < apoptosis_timescale */
    apoptosis_timescale = timescale;
  }
}


/*
 * =============================================================================
 * Start the phagocytosis in the epithelial cell
 * =============================================================================
 */
void
Epithelial::start_phagocytosis (const double time)
{
  is_being_phagocytosed = true;
  phagocytosis_start_time = time;
}


/*
 * =============================================================================
 * End the phagocytosis in the epithelial cell
 * =============================================================================
 */
void
Epithelial::end_phagocytosis ()
{
  is_being_phagocytosed = false;
  phagocytosis_start_time = 0.0;

  change_state(EpithelialState::phagocytosed); 

  /* update the secretions */
  update_secretions();

  flag_bursting = false;
  flag_apoptosis = false;
  flag_phagocytosis = false;
}


/*
 * =============================================================================
 * Set the phagocytosis timescale
 * =============================================================================
 */
void
Epithelial::set_phagocytosis_timescale (const double tms)
{
  (void) tms; // phagocytosis is set by target cell only

  return;

  /* change time scale of apoptosis due to NK action */
  if (CAM_CMP(tms, phagocytosis_timescale, CAM_TOL)) { /* tms < phagocytosis_timescale */
    phagocytosis_timescale = tms;
  }
}
