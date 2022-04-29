#ifndef __EPITHELIAL_H__
#define __EPITHELIAL_H__

#include "../viralrepl/cam_viralrepl.h"
#include "../recepexpr/cam_recepexpr.h"

/* epithelial specific parameters */
struct epithelial_params
{
  /* width and height of epithelial cell */
  double epithelial_width;
  double epithelial_height;

  /*
   * cytokines
   */

  /* T1 IFN */
  double epithelial_healthy_basal_secreted_ifn_1;
  double epithelial_healthy_secreted_ifn_1;
  double epithelial_healthy_uptake_ifn_1;

  double epithelial_infected_basal_secreted_ifn_1;
  double epithelial_infected_secreted_ifn_1;
  double epithelial_infected_uptake_ifn_1;

  double epithelial_ifn_1_EC50;
  double epithelial_ifn_1_hill;
  double epithelial_ifn_1_signal_virus_EC50;
  double epithelial_ifn_1_signal_virus_hill;
  double epithelial_ifn_1_downreg_il_10_EC50;
  double epithelial_ifn_1_downreg_il_10_hill;

  double epithelial_infected_delay_ifn_1;

  /* IL-6 */
  double epithelial_healthy_basal_secreted_il_6;
  double epithelial_healthy_secreted_il_6;
  double epithelial_healthy_uptake_il_6;

  double epithelial_infected_basal_secreted_il_6;
  double epithelial_infected_secreted_il_6;
  double epithelial_infected_uptake_il_6;

  double epithelial_il_6_EC50;
  double epithelial_il_6_hill;
  double epithelial_il_6_signal_virus_EC50;
  double epithelial_il_6_signal_virus_hill;
  double epithelial_il_6_downreg_il_10_EC50;
  double epithelial_il_6_downreg_il_10_hill;

  /* IL-10 */
  double epithelial_healthy_basal_secreted_il_10;
  double epithelial_healthy_secreted_il_10;
  double epithelial_healthy_uptake_il_10;

  double epithelial_infected_basal_secreted_il_10;
  double epithelial_infected_secreted_il_10;
  double epithelial_infected_uptake_il_10;

  double epithelial_il_10_EC50;
  double epithelial_il_10_hill;
  double epithelial_il_10_signal_il_6_EC50;
  double epithelial_il_10_signal_il_6_hill;
  double epithelial_il_10_signal_il_12_EC50;
  double epithelial_il_10_signal_il_12_hill;
  double epithelial_il_10_signal_g_csf_EC50;
  double epithelial_il_10_signal_g_csf_hill;
  double epithelial_il_10_signal_m_csf_EC50;
  double epithelial_il_10_signal_m_csf_hill;

  /* IL-12 */
  double epithelial_healthy_basal_secreted_il_12;
  double epithelial_healthy_secreted_il_12;
  double epithelial_healthy_uptake_il_12;

  double epithelial_infected_basal_secreted_il_12;
  double epithelial_infected_secreted_il_12;
  double epithelial_infected_uptake_il_12;

  double epithelial_il_12_EC50;
  double epithelial_il_12_hill;
  double epithelial_il_12_signal_virus_EC50;
  double epithelial_il_12_signal_virus_hill;
  double epithelial_il_12_downreg_il_10_EC50;
  double epithelial_il_12_downreg_il_10_hill;

  /* G-CSF */
  double epithelial_healthy_basal_secreted_g_csf;
  double epithelial_healthy_secreted_g_csf;
  double epithelial_healthy_uptake_g_csf;

  double epithelial_infected_basal_secreted_g_csf;
  double epithelial_infected_secreted_g_csf;
  double epithelial_infected_uptake_g_csf;

  double epithelial_g_csf_EC50;
  double epithelial_g_csf_hill;
  double epithelial_g_csf_signal_virus_EC50;
  double epithelial_g_csf_signal_virus_hill;
  double epithelial_g_csf_downreg_il_10_EC50;
  double epithelial_g_csf_downreg_il_10_hill;

  /* M-CSF */
  double epithelial_healthy_basal_secreted_m_csf;
  double epithelial_healthy_secreted_m_csf;
  double epithelial_healthy_uptake_m_csf;

  double epithelial_infected_basal_secreted_m_csf;
  double epithelial_infected_secreted_m_csf;
  double epithelial_infected_uptake_m_csf;

  double epithelial_m_csf_EC50;
  double epithelial_m_csf_hill;
  double epithelial_m_csf_signal_virus_EC50;
  double epithelial_m_csf_signal_virus_hill;
  double epithelial_m_csf_downreg_il_10_EC50;
  double epithelial_m_csf_downreg_il_10_hill;

  /* infected epithelial chemokine */
  double epithelial_infected_basal_secreted_epi_inf_ck;
  double epithelial_infected_secreted_epi_inf_ck;
  double epithelial_infected_uptake_epi_inf_ck;

  double epithelial_epi_inf_ck_EC50;
  double epithelial_epi_inf_ck_hill;
  double epithelial_epi_inf_ck_signal_virus_EC50;
  double epithelial_epi_inf_ck_signal_virus_hill;

  /* apoptosis chemokine */
  double epithelial_secreted_apop_ck;

  /* phagocytosis chemokine */
  double epithelial_apoptotic_secreted_phago_ck;

  /* infectious virus - no global memory; defined in replicate */

  /* non-infectious virus - no global memory; defined in replicate */

  /* thresholds */
  double epithelial_burst_threshold_by_virions;

  /* apoptosis */
  double epithelial_apoptosis_timescale;
  double epithelial_apoptosis_max_rate;
  double epithelial_apoptosis_half_max;
  double epithelial_apoptosis_hill;
  double epithelial_apoptosis_by_ifn_max_rate;
  double epithelial_apoptosis_by_ifn_half_max;
  double epithelial_apoptosis_by_ifn_hill;

  /* phagocytosis */
  double epithelial_phagocytosis_timescale;

  /* infectious virus proportion */
  double epithelial_infectious_virus_proportion;

  /* viral replication type */
  const cam_vr_type* vr_typ;

  /* receptor expression type */
  const cam_re_type* re_typ;

};

#endif /* __EPITHELIAL_H__ */
