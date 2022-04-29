#ifndef __MACROPHAGE_H__
#define __MACROPHAGE_H__

/* macrophage specific parameters */
typedef struct {
  /* radius */
  double macrophage_radius;

  /* probability of activating macrophage */
  double macrophage_active_probability;

  /* lifespan */
  double macrophage_resting_lifespan;
  double macrophage_active_lifespan;

  double macrophage_lifespan_downreg_virus_EC50;
  double macrophage_lifespan_downreg_virus_hill;

  /* cytokine secretions */

  /*
   * cytokines
   */

  /* T1 IFN */
  double macrophage_resting_basal_secreted_ifn_1;
  double macrophage_resting_secreted_ifn_1;
  double macrophage_resting_uptake_ifn_1;

  double macrophage_active_basal_secreted_ifn_1;
  double macrophage_active_secreted_ifn_1;
  double macrophage_active_uptake_ifn_1;

  double macrophage_ifn_1_EC50;
  double macrophage_ifn_1_hill;
  double macrophage_ifn_1_signal_virus_EC50;
  double macrophage_ifn_1_signal_virus_hill;
  double macrophage_ifn_1_downreg_il_10_EC50;
  double macrophage_ifn_1_downreg_il_10_hill;

  /* IL-6 */
  double macrophage_resting_basal_secreted_il_6;
  double macrophage_resting_secreted_il_6;
  double macrophage_resting_uptake_il_6;

  double macrophage_active_basal_secreted_il_6;
  double macrophage_active_secreted_il_6;
  double macrophage_active_uptake_il_6;

  double macrophage_il_6_EC50;
  double macrophage_il_6_hill;
  double macrophage_il_6_signal_virus_EC50;
  double macrophage_il_6_signal_virus_hill;
  double macrophage_il_6_downreg_il_10_EC50;
  double macrophage_il_6_downreg_il_10_hill;

  /* IL-10 */
  double macrophage_resting_basal_secreted_il_10;
  double macrophage_resting_secreted_il_10;
  double macrophage_resting_uptake_il_10;

  double macrophage_active_basal_secreted_il_10;
  double macrophage_active_secreted_il_10;
  double macrophage_active_uptake_il_10;

  double macrophage_il_10_EC50;
  double macrophage_il_10_hill;
  double macrophage_il_10_signal_il_6_EC50;
  double macrophage_il_10_signal_il_6_hill;
  double macrophage_il_10_signal_il_12_EC50;
  double macrophage_il_10_signal_il_12_hill;
  double macrophage_il_10_signal_g_csf_EC50;
  double macrophage_il_10_signal_g_csf_hill;
  double macrophage_il_10_signal_m_csf_EC50;
  double macrophage_il_10_signal_m_csf_hill;

  /* IL-12 */
  double macrophage_resting_basal_secreted_il_12;
  double macrophage_resting_secreted_il_12;
  double macrophage_resting_uptake_il_12;

  double macrophage_active_basal_secreted_il_12;
  double macrophage_active_secreted_il_12;
  double macrophage_active_uptake_il_12;

  double macrophage_il_12_EC50;
  double macrophage_il_12_hill;
  double macrophage_il_12_signal_virus_EC50;
  double macrophage_il_12_signal_virus_hill;
  double macrophage_il_12_downreg_il_10_EC50;
  double macrophage_il_12_downreg_il_10_hill;

  /* G-CSF */
  double macrophage_resting_basal_secreted_g_csf;
  double macrophage_resting_secreted_g_csf;
  double macrophage_resting_uptake_g_csf;

  double macrophage_active_basal_secreted_g_csf;
  double macrophage_active_secreted_g_csf;
  double macrophage_active_uptake_g_csf;

  double macrophage_g_csf_EC50;
  double macrophage_g_csf_hill;
  double macrophage_g_csf_signal_virus_EC50;
  double macrophage_g_csf_signal_virus_hill;
  double macrophage_g_csf_downreg_il_10_EC50;
  double macrophage_g_csf_downreg_il_10_hill;

  /* M-CSF */
  double macrophage_resting_basal_secreted_m_csf;
  double macrophage_resting_secreted_m_csf;
  double macrophage_resting_uptake_m_csf;

  double macrophage_active_basal_secreted_m_csf;
  double macrophage_active_secreted_m_csf;
  double macrophage_active_uptake_m_csf;

  double macrophage_m_csf_EC50;
  double macrophage_m_csf_hill;
  double macrophage_m_csf_signal_virus_EC50;
  double macrophage_m_csf_signal_virus_hill;
  double macrophage_m_csf_downreg_il_10_EC50;
  double macrophage_m_csf_downreg_il_10_hill;

  /* infected epithelial chemokine */
  double macrophage_resting_uptake_epi_inf_ck;
  double macrophage_active_uptake_epi_inf_ck;

  /* apoptosis chemokine */
  double macrophage_secreted_apop_ck;

  /* phagocytosis chemokine */
  double macrophage_resting_uptake_phago_ck;
  double macrophage_active_uptake_phago_ck;

  double macrophage_apoptotic_secreted_phago_ck;

  /* infectious virus */
  double macrophage_resting_uptake_virus_inf;
  double macrophage_resting_decay_virus_inf;

  double macrophage_active_uptake_virus_inf;
  double macrophage_active_decay_virus_inf;

  double macrophage_virus_inf_upreg_extracellular_virus_EC50;
  double macrophage_virus_inf_upreg_extracellular_virus_hill;
  double macrophage_virus_inf_downreg_internalised_virus_EC50;
  double macrophage_virus_inf_downreg_internalised_virus_hill;

  /* non-infectious virus - no global memory; not used */

  /* movement */
  double macrophage_resting_movement_rate;
  double macrophage_resting_movement_bias;

  double macrophage_active_movement_rate;
  double macrophage_active_movement_bias;

  double macrophage_movement_upreg_virus_EC50;
  double macrophage_movement_upreg_virus_hill;

  double macrophage_movement_weight_epi_inf_ck;
  double macrophage_movement_weight_phago_ck;

  /* recruitment */
  double macrophage_recruitment_weight_il_6;
  double macrophage_recruitment_weight_il_10;
  double macrophage_recruitment_weight_m_csf;
  double macrophage_recruitment_weight_epi_inf_ck;
  double macrophage_recruitment_weight_phago_ck;

  double macrophage_recruitment_EC50;
  double macrophage_recruitment_hill;

  /* activation and deactivation */
  double macrophage_activate_weight_ifn_1;
  double macrophage_activate_weight_il_6;
  double macrophage_activate_weight_il_10;
  double macrophage_activate_weight_m_csf;
  double macrophage_activate_weight_epi_inf_ck;
  double macrophage_activate_weight_phago_ck;
  double macrophage_activate_weight_virus_inf;

  double macrophage_activate_EC50;
  double macrophage_activate_hill;

  /* bursting */
  double macrophage_resting_burst_threshold;
  double macrophage_active_burst_threshold;

  /* apoptosis */
  double macrophage_apoptosis_timescale;

  /*==========*/
  /* not used */
  double macrophage_apoptosis_max_rate;
  double macrophage_apoptosis_half_max;
  double macrophage_apoptosis_hill;
  /*==========*/

  double macrophage_apoptosis_by_ifn_max_rate;
  double macrophage_apoptosis_by_ifn_half_max;
  double macrophage_apoptosis_by_ifn_hill;


  /* phagocytosis */
  double macrophage_phagocytosis_timescale;
  double macrophage_active_phagocytose_apoptotic_epithelial_probability;
  double macrophage_active_phagocytose_infected_epithelial_probability;
  double macrophage_active_phagocytose_immune_probability;

  bool macrophage_blocking_function;

} macrophage_params;

#endif /* __MACROPHAGE_H__ */
