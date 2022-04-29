#ifndef __MONOCYTE_H__
#define __MONOCYTE_H__

/* monocyte specific parameters */
typedef struct {
  /* radius */
  double monocyte_radius;

  /* probability of activating monocyte */
  double monocyte_active_probability;

  /* lifespan */
  double monocyte_resting_lifespan;
  double monocyte_active_lifespan;

  double monocyte_lifespan_downreg_virus_EC50;
  double monocyte_lifespan_downreg_virus_hill;

  /* cytokine secretions */

  /*
   * cytokines
   */

  /* T1 IFN */
  double monocyte_resting_basal_secreted_ifn_1;
  double monocyte_resting_secreted_ifn_1;
  double monocyte_resting_uptake_ifn_1;

  double monocyte_active_basal_secreted_ifn_1;
  double monocyte_active_secreted_ifn_1;
  double monocyte_active_uptake_ifn_1;

  double monocyte_ifn_1_EC50;
  double monocyte_ifn_1_hill;
  double monocyte_ifn_1_signal_virus_EC50;
  double monocyte_ifn_1_signal_virus_hill;
  double monocyte_ifn_1_downreg_il_10_EC50;
  double monocyte_ifn_1_downreg_il_10_hill;

  /* IL-6 */
  double monocyte_resting_basal_secreted_il_6;
  double monocyte_resting_secreted_il_6;
  double monocyte_resting_uptake_il_6;

  double monocyte_active_basal_secreted_il_6;
  double monocyte_active_secreted_il_6;
  double monocyte_active_uptake_il_6;

  double monocyte_il_6_EC50;
  double monocyte_il_6_hill;
  double monocyte_il_6_signal_virus_EC50;
  double monocyte_il_6_signal_virus_hill;
  double monocyte_il_6_downreg_il_10_EC50;
  double monocyte_il_6_downreg_il_10_hill;

  /* IL-10 */
  double monocyte_resting_basal_secreted_il_10;
  double monocyte_resting_secreted_il_10;
  double monocyte_resting_uptake_il_10;

  double monocyte_active_basal_secreted_il_10;
  double monocyte_active_secreted_il_10;
  double monocyte_active_uptake_il_10;

  double monocyte_il_10_EC50;
  double monocyte_il_10_hill;
  double monocyte_il_10_signal_il_6_EC50;
  double monocyte_il_10_signal_il_6_hill;
  double monocyte_il_10_signal_il_12_EC50;
  double monocyte_il_10_signal_il_12_hill;
  double monocyte_il_10_signal_g_csf_EC50;
  double monocyte_il_10_signal_g_csf_hill;
  double monocyte_il_10_signal_m_csf_EC50;
  double monocyte_il_10_signal_m_csf_hill;

  /* IL-12 */
  double monocyte_resting_basal_secreted_il_12;
  double monocyte_resting_secreted_il_12;
  double monocyte_resting_uptake_il_12;

  double monocyte_active_basal_secreted_il_12;
  double monocyte_active_secreted_il_12;
  double monocyte_active_uptake_il_12;

  double monocyte_il_12_EC50;
  double monocyte_il_12_hill;
  double monocyte_il_12_signal_virus_EC50;
  double monocyte_il_12_signal_virus_hill;
  double monocyte_il_12_downreg_il_10_EC50;
  double monocyte_il_12_downreg_il_10_hill;

  /* G-CSF */
  double monocyte_resting_basal_secreted_g_csf;
  double monocyte_resting_secreted_g_csf;
  double monocyte_resting_uptake_g_csf;

  double monocyte_active_basal_secreted_g_csf;
  double monocyte_active_secreted_g_csf;
  double monocyte_active_uptake_g_csf;

  double monocyte_g_csf_EC50;
  double monocyte_g_csf_hill;
  double monocyte_g_csf_signal_virus_EC50;
  double monocyte_g_csf_signal_virus_hill;
  double monocyte_g_csf_downreg_il_10_EC50;
  double monocyte_g_csf_downreg_il_10_hill;

  /* M-CSF */
  double monocyte_resting_basal_secreted_m_csf;
  double monocyte_resting_secreted_m_csf;
  double monocyte_resting_uptake_m_csf;

  double monocyte_active_basal_secreted_m_csf;
  double monocyte_active_secreted_m_csf;
  double monocyte_active_uptake_m_csf;

  double monocyte_m_csf_EC50;
  double monocyte_m_csf_hill;
  double monocyte_m_csf_signal_virus_EC50;
  double monocyte_m_csf_signal_virus_hill;
  double monocyte_m_csf_downreg_il_10_EC50;
  double monocyte_m_csf_downreg_il_10_hill;

  /* infected epithelial chemokine */
  double monocyte_resting_uptake_epi_inf_ck;
  double monocyte_active_uptake_epi_inf_ck;

  /* apoptosis chemokine */
  double monocyte_secreted_apop_ck;

  /* phagocytosis chemokine */
  double monocyte_resting_uptake_phago_ck;
  double monocyte_active_uptake_phago_ck;

  double monocyte_apoptotic_secreted_phago_ck;

  /* infectious virus */
  double monocyte_resting_uptake_virus_inf;
  double monocyte_resting_decay_virus_inf;

  double monocyte_active_uptake_virus_inf;
  double monocyte_active_decay_virus_inf;

  double monocyte_virus_inf_upreg_extracellular_virus_EC50;
  double monocyte_virus_inf_upreg_extracellular_virus_hill;
  double monocyte_virus_inf_downreg_internalised_virus_EC50;
  double monocyte_virus_inf_downreg_internalised_virus_hill;

  /* non-infectious virus - no global memory; not used */

  /* movement */
  double monocyte_resting_movement_rate;
  double monocyte_resting_movement_bias;

  double monocyte_active_movement_rate;
  double monocyte_active_movement_bias;

  double monocyte_movement_upreg_virus_EC50;
  double monocyte_movement_upreg_virus_hill;

  double monocyte_movement_weight_epi_inf_ck;
  double monocyte_movement_weight_phago_ck;

  /* recruitment */
  double monocyte_recruitment_weight_il_6;
  double monocyte_recruitment_weight_il_10;
  double monocyte_recruitment_weight_m_csf;
  double monocyte_recruitment_weight_epi_inf_ck;
  double monocyte_recruitment_weight_phago_ck;

  double monocyte_recruitment_EC50;
  double monocyte_recruitment_hill;

  /* activation and deactivation */
  double monocyte_activate_weight_ifn_1;
  double monocyte_activate_weight_il_6;
  double monocyte_activate_weight_il_10;
  double monocyte_activate_weight_m_csf;
  double monocyte_activate_weight_epi_inf_ck;
  double monocyte_activate_weight_phago_ck;
  double monocyte_activate_weight_virus_inf;

  double monocyte_activate_EC50;
  double monocyte_activate_hill;

  /* bursting */
  double monocyte_resting_burst_threshold;
  double monocyte_active_burst_threshold;

  /* apoptosis */
  double monocyte_apoptosis_timescale;
  double monocyte_apoptosis_max_rate;
  double monocyte_apoptosis_half_max;
  double monocyte_apoptosis_hill;
  double monocyte_apoptosis_by_ifn_max_rate;
  double monocyte_apoptosis_by_ifn_half_max;
  double monocyte_apoptosis_by_ifn_hill;

  /* phagocytosis */
  double monocyte_phagocytosis_timescale;
  double monocyte_active_phagocytose_apoptotic_epithelial_probability;
  double monocyte_active_phagocytose_infected_epithelial_probability;
  double monocyte_active_phagocytose_immune_probability;

  /* conversion */
  double monocyte_active_to_mac_probability;

  bool monocyte_blocking_function;

} monocyte_params;

#endif /* __MONOCYTE_H__ */
