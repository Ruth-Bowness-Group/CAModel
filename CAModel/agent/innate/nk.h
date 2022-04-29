#ifndef __NK_H__
#define __NK_H__

/* nk specific parameters */
typedef struct {
  /* radius */
  double nk_radius;

  /* probability of activating nk */
  double nk_active_probability;

  /* lifespan */
  double nk_resting_lifespan;
  double nk_active_lifespan;

  double nk_lifespan_downreg_virus_EC50;
  double nk_lifespan_downreg_virus_hill;

  /* cytokine secretions */

  /*
   * cytokines
   */

  /* T1 IFN */
  double nk_resting_basal_secreted_ifn_1;
  double nk_resting_secreted_ifn_1;
  double nk_resting_uptake_ifn_1;

  double nk_active_basal_secreted_ifn_1;
  double nk_active_secreted_ifn_1;
  double nk_active_uptake_ifn_1;

  double nk_ifn_1_EC50;
  double nk_ifn_1_hill;
  double nk_ifn_1_signal_virus_EC50;
  double nk_ifn_1_signal_virus_hill;
  double nk_ifn_1_downreg_il_10_EC50;
  double nk_ifn_1_downreg_il_10_hill;

  /* IL-6 */
  double nk_resting_basal_secreted_il_6;
  double nk_resting_secreted_il_6;
  double nk_resting_uptake_il_6;

  double nk_active_basal_secreted_il_6;
  double nk_active_secreted_il_6;
  double nk_active_uptake_il_6;

  double nk_il_6_EC50;
  double nk_il_6_hill;
  double nk_il_6_signal_virus_EC50;
  double nk_il_6_signal_virus_hill;
  double nk_il_6_downreg_il_10_EC50;
  double nk_il_6_downreg_il_10_hill;

  /* IL-10 */
  double nk_resting_basal_secreted_il_10;
  double nk_resting_secreted_il_10;
  double nk_resting_uptake_il_10;

  double nk_active_basal_secreted_il_10;
  double nk_active_secreted_il_10;
  double nk_active_uptake_il_10;

  double nk_il_10_EC50;
  double nk_il_10_hill;
  double nk_il_10_signal_il_6_EC50;
  double nk_il_10_signal_il_6_hill;
  double nk_il_10_signal_il_12_EC50;
  double nk_il_10_signal_il_12_hill;
  double nk_il_10_signal_g_csf_EC50;
  double nk_il_10_signal_g_csf_hill;
  double nk_il_10_signal_m_csf_EC50;
  double nk_il_10_signal_m_csf_hill;

  /* IL-12 */
  double nk_resting_basal_secreted_il_12;
  double nk_resting_secreted_il_12;
  double nk_resting_uptake_il_12;

  double nk_active_basal_secreted_il_12;
  double nk_active_secreted_il_12;
  double nk_active_uptake_il_12;

  double nk_il_12_EC50;
  double nk_il_12_hill;
  double nk_il_12_signal_virus_EC50;
  double nk_il_12_signal_virus_hill;
  double nk_il_12_downreg_il_10_EC50;
  double nk_il_12_downreg_il_10_hill;

  /* G-CSF */
  double nk_resting_basal_secreted_g_csf;
  double nk_resting_secreted_g_csf;
  double nk_resting_uptake_g_csf;

  double nk_active_basal_secreted_g_csf;
  double nk_active_secreted_g_csf;
  double nk_active_uptake_g_csf;

  double nk_g_csf_EC50;
  double nk_g_csf_hill;
  double nk_g_csf_signal_virus_EC50;
  double nk_g_csf_signal_virus_hill;
  double nk_g_csf_downreg_il_10_EC50;
  double nk_g_csf_downreg_il_10_hill;

  /* M-CSF */
  double nk_resting_basal_secreted_m_csf;
  double nk_resting_secreted_m_csf;
  double nk_resting_uptake_m_csf;

  double nk_active_basal_secreted_m_csf;
  double nk_active_secreted_m_csf;
  double nk_active_uptake_m_csf;

  double nk_m_csf_EC50;
  double nk_m_csf_hill;
  double nk_m_csf_signal_virus_EC50;
  double nk_m_csf_signal_virus_hill;
  double nk_m_csf_downreg_il_10_EC50;
  double nk_m_csf_downreg_il_10_hill;

  /* infected epithelial chemokine */
  double nk_resting_uptake_epi_inf_ck;
  double nk_active_uptake_epi_inf_ck;

  /* apoptosis chemokine */
  double nk_secreted_apop_ck;

  double nk_resting_uptake_apop_ck;
  double nk_active_uptake_apop_ck;

  /* phagocytosis chemokine */
  double nk_apoptotic_secreted_phago_ck;

  /* infectious virus - no global memory; not used */

  /* non-infectious virus - no global memory; not used */

  /* movement */
  double nk_resting_movement_rate;
  double nk_resting_movement_bias;

  double nk_active_movement_rate;
  double nk_active_movement_bias;

  double nk_movement_upreg_virus_EC50;
  double nk_movement_upreg_virus_hill;

  double nk_movement_weight_epi_inf_ck;
  double nk_movement_weight_apop_ck;

  /* recruitment */
  double nk_recruitment_weight_il_6;
  double nk_recruitment_weight_il_10;
  double nk_recruitment_weight_il_12;
  double nk_recruitment_weight_epi_inf_ck;
  double nk_recruitment_weight_apop_ck;

  double nk_recruitment_EC50;
  double nk_recruitment_hill;

  /* activation and deactivation */
  double nk_activate_weight_ifn_1;
  double nk_activate_weight_il_6;
  double nk_activate_weight_il_10;
  double nk_activate_weight_il_12;
  double nk_activate_weight_epi_inf_ck;
  double nk_activate_weight_apop_ck;
  double nk_activate_weight_virus_inf;

  double nk_activate_EC50;
  double nk_activate_hill;

  /* apoptosis */
  double nk_apoptosis_timescale;
  double nk_apoptosis_max_rate;
  double nk_apoptosis_half_max;
  double nk_apoptosis_hill;
  double nk_apoptosis_by_ifn_max_rate;
  double nk_apoptosis_by_ifn_half_max;
  double nk_apoptosis_by_ifn_hill;
  double nk_active_apoptose_epithelial_probability;
  double nk_active_apoptose_immune_probability;

  /* phagocytosis */
  double nk_phagocytosis_timescale;

} nk_params;

#endif /* __NK_H__ */
