#ifndef __REACTDIFF_H__
#define __REACTDIFF_H__


/* struct containing all of the reaction diffusion specific parameters */
struct rd_params
{
  /* spatial step of the grid */
  double rd_spatial_step;

  /* T1 IFN */
  double rd_ifn_1_diffusion;
  double rd_ifn_1_decay;

  /* IL-6 */
  double rd_il_6_diffusion;
  double rd_il_6_decay;

  /* IL-10 */
  double rd_il_10_diffusion;
  double rd_il_10_decay;

  /* IL-12 */
  double rd_il_12_diffusion;
  double rd_il_12_decay;

  /* G-CSF */
  double rd_g_csf_diffusion;
  double rd_g_csf_decay;

  /* M-CSF */
  double rd_m_csf_diffusion;
  double rd_m_csf_decay;

  /* infected epithelial chemokine */
  double rd_epi_inf_ck_diffusion;
  double rd_epi_inf_ck_decay;

  /* apoptosis chemokine */
  double rd_apop_ck_diffusion;
  double rd_apop_ck_decay;

  /* phagocytosis chemokine */
  double rd_phago_ck_diffusion;
  double rd_phago_ck_decay;

  /* infectious virus */
  double rd_virus_inf_diffusion;
  double rd_virus_inf_decay;

  /* non-infectious virus */
  double rd_virus_ninf_diffusion;
  double rd_virus_ninf_decay;

};

#endif /* __REACTDIFF_H__ */
