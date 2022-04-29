#ifndef __NEUTROPHIL_H__
#define __NEUTROPHIL_H__

/* neutrophil specific parameters */
typedef struct {
  /* radius */
  double neutrophil_radius;

  /* probability of activating neutrophil */
  double neutrophil_active_probability;

  /* lifespan */
  double neutrophil_resting_lifespan;
  double neutrophil_active_lifespan;

  double neutrophil_lifespan_downreg_virus_EC50;
  double neutrophil_lifespan_downreg_virus_hill;

  /* cytokine secretions */

  /*
   * cytokines
   */

  /* T1 IFN */
  double neutrophil_resting_basal_secreted_ifn_1;
  double neutrophil_resting_secreted_ifn_1;
  double neutrophil_resting_uptake_ifn_1;

  double neutrophil_active_basal_secreted_ifn_1;
  double neutrophil_active_secreted_ifn_1;
  double neutrophil_active_uptake_ifn_1;

  double neutrophil_ifn_1_EC50;
  double neutrophil_ifn_1_hill;
  double neutrophil_ifn_1_signal_virus_EC50;
  double neutrophil_ifn_1_signal_virus_hill;
  double neutrophil_ifn_1_downreg_il_10_EC50;
  double neutrophil_ifn_1_downreg_il_10_hill;

  /* IL-6 */
  double neutrophil_resting_basal_secreted_il_6;
  double neutrophil_resting_secreted_il_6;
  double neutrophil_resting_uptake_il_6;

  double neutrophil_active_basal_secreted_il_6;
  double neutrophil_active_secreted_il_6;
  double neutrophil_active_uptake_il_6;

  double neutrophil_il_6_EC50;
  double neutrophil_il_6_hill;
  double neutrophil_il_6_signal_virus_EC50;
  double neutrophil_il_6_signal_virus_hill;
  double neutrophil_il_6_downreg_il_10_EC50;
  double neutrophil_il_6_downreg_il_10_hill;

  /* IL-10 */
  double neutrophil_resting_basal_secreted_il_10;
  double neutrophil_resting_secreted_il_10;
  double neutrophil_resting_uptake_il_10;

  double neutrophil_active_basal_secreted_il_10;
  double neutrophil_active_secreted_il_10;
  double neutrophil_active_uptake_il_10;

  double neutrophil_il_10_EC50;
  double neutrophil_il_10_hill;
  double neutrophil_il_10_signal_il_6_EC50;
  double neutrophil_il_10_signal_il_6_hill;
  double neutrophil_il_10_signal_il_12_EC50;
  double neutrophil_il_10_signal_il_12_hill;
  double neutrophil_il_10_signal_g_csf_EC50;
  double neutrophil_il_10_signal_g_csf_hill;
  double neutrophil_il_10_signal_m_csf_EC50;
  double neutrophil_il_10_signal_m_csf_hill;

  /* IL-12 */
  double neutrophil_resting_basal_secreted_il_12;
  double neutrophil_resting_secreted_il_12;
  double neutrophil_resting_uptake_il_12;

  double neutrophil_active_basal_secreted_il_12;
  double neutrophil_active_secreted_il_12;
  double neutrophil_active_uptake_il_12;

  double neutrophil_il_12_EC50;
  double neutrophil_il_12_hill;
  double neutrophil_il_12_signal_virus_EC50;
  double neutrophil_il_12_signal_virus_hill;
  double neutrophil_il_12_downreg_il_10_EC50;
  double neutrophil_il_12_downreg_il_10_hill;

  /* G-CSF */
  double neutrophil_resting_basal_secreted_g_csf;
  double neutrophil_resting_secreted_g_csf;
  double neutrophil_resting_uptake_g_csf;

  double neutrophil_active_basal_secreted_g_csf;
  double neutrophil_active_secreted_g_csf;
  double neutrophil_active_uptake_g_csf;

  double neutrophil_g_csf_EC50;
  double neutrophil_g_csf_hill;
  double neutrophil_g_csf_signal_virus_EC50;
  double neutrophil_g_csf_signal_virus_hill;
  double neutrophil_g_csf_downreg_il_10_EC50;
  double neutrophil_g_csf_downreg_il_10_hill;

  /* M-CSF */
  double neutrophil_resting_basal_secreted_m_csf;
  double neutrophil_resting_secreted_m_csf;
  double neutrophil_resting_uptake_m_csf;

  double neutrophil_active_basal_secreted_m_csf;
  double neutrophil_active_secreted_m_csf;
  double neutrophil_active_uptake_m_csf;

  double neutrophil_m_csf_EC50;
  double neutrophil_m_csf_hill;
  double neutrophil_m_csf_signal_virus_EC50;
  double neutrophil_m_csf_signal_virus_hill;
  double neutrophil_m_csf_downreg_il_10_EC50;
  double neutrophil_m_csf_downreg_il_10_hill;

  /* infected epithelial chemokine */
  double neutrophil_resting_uptake_epi_inf_ck;
  double neutrophil_active_uptake_epi_inf_ck;

  /* apoptosis chemokine */
  double neutrophil_secreted_apop_ck;

  /* phagocytosis chemokine */
  double neutrophil_resting_uptake_phago_ck;
  double neutrophil_active_uptake_phago_ck;

  double neutrophil_apoptotic_secreted_phago_ck;

  /* infectious virus */
  double neutrophil_resting_uptake_virus_inf;
  double neutrophil_resting_decay_virus_inf;

  double neutrophil_active_uptake_virus_inf;
  double neutrophil_active_decay_virus_inf;

  double neutrophil_virus_inf_upreg_extracellular_virus_EC50;
  double neutrophil_virus_inf_upreg_extracellular_virus_hill;
  double neutrophil_virus_inf_downreg_internalised_virus_EC50;
  double neutrophil_virus_inf_downreg_internalised_virus_hill;

  /* non-infectious virus - no global memory; not used */

  /* movement */
  double neutrophil_resting_movement_rate;
  double neutrophil_resting_movement_bias;

  double neutrophil_active_movement_rate;
  double neutrophil_active_movement_bias;

  double neutrophil_movement_upreg_virus_EC50;
  double neutrophil_movement_upreg_virus_hill;

  double neutrophil_movement_weight_epi_inf_ck;
  double neutrophil_movement_weight_phago_ck;

  /* recruitment */
  double neutrophil_recruitment_weight_il_6;
  double neutrophil_recruitment_weight_il_10;
  double neutrophil_recruitment_weight_g_csf;
  double neutrophil_recruitment_weight_epi_inf_ck;
  double neutrophil_recruitment_weight_phago_ck;

  double neutrophil_recruitment_EC50;
  double neutrophil_recruitment_hill;

  /* activation and deactivation */
  double neutrophil_activate_weight_ifn_1;
  double neutrophil_activate_weight_il_6;
  double neutrophil_activate_weight_il_10;
  double neutrophil_activate_weight_g_csf;
  double neutrophil_activate_weight_epi_inf_ck;
  double neutrophil_activate_weight_phago_ck;
  double neutrophil_activate_weight_virus_inf;

  double neutrophil_activate_EC50;
  double neutrophil_activate_hill;

  /* bursting */
  double neutrophil_resting_burst_threshold;
  double neutrophil_active_burst_threshold;

  /* apoptosis */
  double neutrophil_apoptosis_timescale;
  double neutrophil_apoptosis_max_rate;
  double neutrophil_apoptosis_half_max;
  double neutrophil_apoptosis_hill;
  double neutrophil_apoptosis_by_ifn_max_rate;
  double neutrophil_apoptosis_by_ifn_half_max;
  double neutrophil_apoptosis_by_ifn_hill;

  /* phagocytosis */
  double neutrophil_phagocytosis_timescale;
  double neutrophil_active_phagocytose_apoptotic_epithelial_probability;
  double neutrophil_active_phagocytose_infected_epithelial_probability;
  double neutrophil_active_phagocytose_immune_probability;

  bool neutrophil_blocking_function;

} neutrophil_params;

#endif /* __NEUTROPHIL_H__ */
