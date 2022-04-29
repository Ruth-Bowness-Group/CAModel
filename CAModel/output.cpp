#include "cam_output.h"
#include "cam_error.h"


/*
 * =============================================================================
 * Open the output files
 * =============================================================================
 */
int
cam_output_open (const std::string output_folder, struct cam_output_files* files)
{
  std::string outf;

  /* global constants */
  outf = output_folder + std::string("/consts.txt");
  files->consts = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->consts, "consts.txt");

  /*
   * ---------------------------------------------------------------------------
   * grid output codes
   * ---------------------------------------------------------------------------
   */

  /* location grid codes */
  outf = output_folder + std::string("/grid_codes.txt");
  files->grid_codes = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->grid_codes, "grid_codes.txt");


  /* macrophage grid codes */
  outf = output_folder + std::string("/grid_codes_macrophage_alveolar_resting.txt");
  files->grid_codes_macrophage_alveolar_resting = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->grid_codes_macrophage_alveolar_resting, "grid_codes_macrophage_alveolar_resting.txt");

  outf = output_folder + std::string("/grid_codes_macrophage_alveolar_active.txt");
  files->grid_codes_macrophage_alveolar_active = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->grid_codes_macrophage_alveolar_active, "grid_codes_macrophage_alveolar_active.txt");

  outf = output_folder + std::string("/grid_codes_macrophage_alveolar_apoptotic.txt");
  files->grid_codes_macrophage_alveolar_apoptotic = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->grid_codes_macrophage_alveolar_apoptotic, "grid_codes_macrophage_alveolar_apoptotic.txt");


  /* nk grid codes */
  outf = output_folder + std::string("/grid_codes_nk_cytotoxic_resting.txt");
  files->grid_codes_nk_cytotoxic_resting = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->grid_codes_nk_cytotoxic_resting, "grid_codes_nk_cytotoxic_resting.txt");

  outf = output_folder + std::string("/grid_codes_nk_cytotoxic_active.txt");
  files->grid_codes_nk_cytotoxic_active = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->grid_codes_nk_cytotoxic_active, "grid_codes_nk_cytotoxic_active.txt");

  outf = output_folder + std::string("/grid_codes_nk_cytotoxic_apoptotic.txt");
  files->grid_codes_nk_cytotoxic_apoptotic = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->grid_codes_nk_cytotoxic_apoptotic, "grid_codes_nk_cytotoxic_apoptotic.txt");


  /* neutrophil grid codes */
  outf = output_folder + std::string("/grid_codes_neutrophil_null_resting.txt");
  files->grid_codes_neutrophil_null_resting = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->grid_codes_neutrophil_null_resting, "grid_codes_neutrophil_null_resting.txt");

  outf = output_folder + std::string("/grid_codes_neutrophil_null_active.txt");
  files->grid_codes_neutrophil_null_active = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->grid_codes_neutrophil_null_active, "grid_codes_neutrophil_null_active.txt");

  outf = output_folder + std::string("/grid_codes_neutrophil_null_apoptotic.txt");
  files->grid_codes_neutrophil_null_apoptotic = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->grid_codes_neutrophil_null_apoptotic, "grid_codes_neutrophil_null_apoptotic.txt");


  /* monocyte grid codes */
  outf = output_folder + std::string("/grid_codes_monocyte_null_resting.txt");
  files->grid_codes_monocyte_null_resting = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->grid_codes_monocyte_null_resting, "grid_codes_monocyte_null_resting.txt");

  outf = output_folder + std::string("/grid_codes_monocyte_null_active.txt");
  files->grid_codes_monocyte_null_active = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->grid_codes_monocyte_null_active, "grid_codes_monocyte_null_active.txt");

  outf = output_folder + std::string("/grid_codes_monocyte_null_apoptotic.txt");
  files->grid_codes_monocyte_null_apoptotic = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->grid_codes_monocyte_null_apoptotic, "grid_codes_monocyte_null_apoptotic.txt");


  /*
   * ---------------------------------------------------------------------------
   * Cell counts
   * ---------------------------------------------------------------------------
   */

  /* epithelial cell counts */
  outf = output_folder + std::string("/cell_counts_epithelial_healthy.txt");
  files->cell_counts_epithelial_healthy = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_epithelial_healthy, "cell_counts_epithelial_healthy.txt");

  outf = output_folder + std::string("/cell_counts_epithelial_infected.txt");
  files->cell_counts_epithelial_infected = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_epithelial_infected, "cell_counts_epithelial_infected.txt");

  outf = output_folder + std::string("/cell_counts_epithelial_apoptotic.txt");
  files->cell_counts_epithelial_apoptotic = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_epithelial_apoptotic, "cell_counts_epithelial_apoptotic.txt");

  outf = output_folder + std::string("/cell_counts_epithelial_phagocytosed.txt");
  files->cell_counts_epithelial_phagocytosed = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_epithelial_phagocytosed, "cell_counts_epithelial_phagocytosed.txt");

  outf = output_folder + std::string("/cell_counts_epithelial_burst.txt");
  files->cell_counts_epithelial_burst = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_epithelial_burst, "cell_counts_epithelial_burst.txt");

  outf = output_folder + std::string("/cell_counts_epithelial_ratio.txt");
  files->cell_counts_epithelial_ratio = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_epithelial_ratio, "cell_counts_epithelial_ratio.txt");


  /* infectious virus cell counts */
  outf = output_folder + std::string("/cell_counts_virus_inf_extracellular.txt");
  files->cell_counts_virus_inf_extracellular = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_virus_inf_extracellular, "cell_counts_virus_inf_extracellular.txt");

  outf = output_folder + std::string("/cell_counts_virus_inf_intracellular.txt");
  files->cell_counts_virus_inf_intracellular = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_virus_inf_intracellular, "cell_counts_virus_inf_intracellular.txt");

  outf = output_folder + std::string("/cell_counts_virus_inf_total.txt");
  files->cell_counts_virus_inf_total = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_virus_inf_total, "cell_counts_virus_inf_total.txt");


  /* non-infectious virus cell counts */
  outf = output_folder + std::string("/cell_counts_virus_ninf_extracellular.txt");
  files->cell_counts_virus_ninf_extracellular = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_virus_ninf_extracellular, "cell_counts_virus_ninf_extracellular.txt");

  outf = output_folder + std::string("/cell_counts_virus_ninf_total.txt");
  files->cell_counts_virus_ninf_total = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_virus_ninf_total, "cell_counts_virus_ninf_total.txt");


  /* virus cell counts */
  outf = output_folder + std::string("/cell_counts_virus_extracellular_total.txt");
  files->cell_counts_virus_extracellular_total = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_virus_extracellular_total, "cell_counts_virus_extracellular_total.txt");

  outf = output_folder + std::string("/cell_counts_virus_total.txt");
  files->cell_counts_virus_total = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_virus_total, "cell_counts_virus_total.txt");


  /* macrophage cell counts */
  outf = output_folder + std::string("/cell_counts_macrophages_alveolar_resting.txt");
  files->cell_counts_macrophages_alveolar_resting = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_macrophages_alveolar_resting, "cell_counts_macrophages_alveolar_resting.txt");

  outf = output_folder + std::string("/cell_counts_macrophages_alveolar_active.txt");
  files->cell_counts_macrophages_alveolar_active = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_macrophages_alveolar_active, "cell_counts_macrophages_alveolar_active.txt");

  outf = output_folder + std::string("/cell_counts_macrophages_alveolar_apoptotic.txt");
  files->cell_counts_macrophages_alveolar_apoptotic = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_macrophages_alveolar_apoptotic, "cell_counts_macrophages_alveolar_apoptotic.txt");

  outf = output_folder + std::string("/cell_counts_macrophages_alveolar.txt");
  files->cell_counts_macrophages_alveolar = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_macrophages_alveolar, "cell_counts_macrophages_alveolar.txt");


  /* nk cell counts */
  outf = output_folder + std::string("/cell_counts_nks_cytotoxic_resting.txt");
  files->cell_counts_nks_cytotoxic_resting = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_nks_cytotoxic_resting, "cell_counts_nks_cytotoxic_resting.txt");

  outf = output_folder + std::string("/cell_counts_nks_cytotoxic_active.txt");
  files->cell_counts_nks_cytotoxic_active = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_nks_cytotoxic_active, "cell_counts_nks_cytotoxic_active.txt");

  outf = output_folder + std::string("/cell_counts_nks_cytotoxic_apoptotic.txt");
  files->cell_counts_nks_cytotoxic_apoptotic = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_nks_cytotoxic_apoptotic, "cell_counts_nks_cytotoxic_apoptotic.txt");

  outf = output_folder + std::string("/cell_counts_nks_cytotoxic.txt");
  files->cell_counts_nks_cytotoxic = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_nks_cytotoxic, "cell_counts_nks_cytotoxic.txt");


  /* neutrophil cell counts */
  outf = output_folder + std::string("/cell_counts_neutrophils_null_resting.txt");
  files->cell_counts_neutrophils_null_resting = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_neutrophils_null_resting, "cell_counts_neutrophils_null_resting.txt");

  outf = output_folder + std::string("/cell_counts_neutrophils_null_active.txt");
  files->cell_counts_neutrophils_null_active = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_neutrophils_null_active, "cell_counts_neutrophils_null_active.txt");

  outf = output_folder + std::string("/cell_counts_neutrophils_null_apoptotic.txt");
  files->cell_counts_neutrophils_null_apoptotic = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_neutrophils_null_apoptotic, "cell_counts_neutrophils_null_apoptotic.txt");

  outf = output_folder + std::string("/cell_counts_neutrophils_null.txt");
  files->cell_counts_neutrophils_null = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_neutrophils_null, "cell_counts_neutrophils_null.txt");


  /* monocyte cell counts */
  outf = output_folder + std::string("/cell_counts_monocytes_null_resting.txt");
  files->cell_counts_monocytes_null_resting = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_monocytes_null_resting, "cell_counts_monocytes_null_resting.txt");

  outf = output_folder + std::string("/cell_counts_monocytes_null_active.txt");
  files->cell_counts_monocytes_null_active = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_monocytes_null_active, "cell_counts_monocytes_null_active.txt");

  outf = output_folder + std::string("/cell_counts_monocytes_null_apoptotic.txt");
  files->cell_counts_monocytes_null_apoptotic = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_monocytes_null_apoptotic, "cell_counts_monocytes_null_apoptotic.txt");

  outf = output_folder + std::string("/cell_counts_monocytes_null.txt");
  files->cell_counts_monocytes_null = fopen(outf.c_str(), "w");
  CAM_ERROR_CHECK_IO(files->cell_counts_monocytes_null, "cell_counts_monocytes_null.txt");


  /*
   * ---------------------------------------------------------------------------
   * reaction diffusion systems
   * ---------------------------------------------------------------------------
   */

  /* type I interferon */
  outf = output_folder + std::string("/rd_ifn_1.txt");
  files->rd_ifn_1 = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_ifn_1, "rd_ifn_1.txt");

  outf = output_folder + std::string("/rd_total_ifn_1.txt");
  files->rd_total_ifn_1 = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_total_ifn_1, "rd_total_ifn_1.txt");

  /* il-6 */
  outf = output_folder + std::string("/rd_il_6.txt");
  files->rd_il_6 = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_il_6, "rd_il_6.txt");

  outf = output_folder + std::string("/rd_total_il_6.txt");
  files->rd_total_il_6 = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_total_il_6, "rd_total_il_6.txt");

  /* il-10 */
  outf = output_folder + std::string("/rd_il_10.txt");
  files->rd_il_10 = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_il_10, "rd_il_10.txt");

  outf = output_folder + std::string("/rd_total_il_10.txt");
  files->rd_total_il_10 = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_total_il_10, "rd_total_il_10.txt");

  /* il-12 */
  outf = output_folder + std::string("/rd_il_12.txt");
  files->rd_il_12 = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_il_12, "rd_il_12.txt");

  outf = output_folder + std::string("/rd_total_il_12.txt");
  files->rd_total_il_12 = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_total_il_12, "rd_total_il_12.txt");

  /* g_csf */
  outf = output_folder + std::string("/rd_g_csf.txt");
  files->rd_g_csf = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_g_csf, "rd_g_csf.txt");

  outf = output_folder + std::string("/rd_total_g_csf.txt");
  files->rd_total_g_csf = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_total_g_csf, "rd_total_g_csf.txt");

  /* m_csf */
  outf = output_folder + std::string("/rd_m_csf.txt");
  files->rd_m_csf = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_m_csf, "rd_m_csf.txt");

  outf = output_folder + std::string("/rd_total_m_csf.txt");
  files->rd_total_m_csf = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_total_m_csf, "rd_total_m_csf.txt");

  /* epi_inf_ck */
  outf = output_folder + std::string("/rd_epi_inf_ck.txt");
  files->rd_epi_inf_ck = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_epi_inf_ck, "rd_epi_inf_ck.txt");

  outf = output_folder + std::string("/rd_total_epi_inf_ck.txt");
  files->rd_total_epi_inf_ck = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_total_epi_inf_ck, "rd_total_epi_inf_ck.txt");

  /* apoptosis chemokine */
  outf = output_folder + std::string("/rd_apop_ck.txt");
  files->rd_apop_ck = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_apop_ck, "rd_apop_ck.txt");

  outf = output_folder + std::string("/rd_total_apop_ck.txt");
  files->rd_total_apop_ck = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_total_apop_ck, "rd_total_apop_ck.txt");

  /* phagocytosis chemokine */
  outf = output_folder + std::string("/rd_phago_ck.txt");
  files->rd_phago_ck = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_phago_ck, "rd_phago_ck.txt");

  outf = output_folder + std::string("/rd_total_phago_ck.txt");
  files->rd_total_phago_ck = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_total_phago_ck, "rd_total_phago_ck.txt");

  /* infectious virus */
  outf = output_folder + std::string("/rd_virus_inf.txt");
  files->rd_virus_inf = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_virus_inf, "rd_virus_inf.txt");

  outf = output_folder + std::string("/rd_total_virus_inf.txt");
  files->rd_total_virus_inf = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_total_virus_inf, "rd_total_virus_inf.txt");

  /* non-infectious virus */
  outf = output_folder + std::string("/rd_virus_ninf.txt");
  files->rd_virus_ninf = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_virus_ninf, "rd_virus_ninf.txt");

  outf = output_folder + std::string("/rd_total_virus_ninf.txt");
  files->rd_total_virus_ninf = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_total_virus_ninf, "rd_total_virus_ninf.txt");

  /* virus */
  outf = output_folder + std::string("/rd_total_virus.txt");
  files->rd_total_virus = fopen(outf.c_str(),"w");
  CAM_ERROR_CHECK_IO(files->rd_total_virus, "rd_total_virus.txt");

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Close the output files
 * =============================================================================
 */
void
cam_output_close (struct cam_output_files* files)
{
  /* constants */
  fclose(files->consts);

  /*
   * ---------------------------------------------------------------------------
   * grid output codes
   * ---------------------------------------------------------------------------
   */

  /* location grid codes */
  fclose(files->grid_codes);

  /* macrophage grid codes */
  fclose(files->grid_codes_macrophage_alveolar_resting);
  fclose(files->grid_codes_macrophage_alveolar_active);
  fclose(files->grid_codes_macrophage_alveolar_apoptotic);

  /* nk grid codes */
  fclose(files->grid_codes_nk_cytotoxic_resting);
  fclose(files->grid_codes_nk_cytotoxic_active);
  fclose(files->grid_codes_nk_cytotoxic_apoptotic);

  /* neutrophil grid codes */
  fclose(files->grid_codes_neutrophil_null_resting);
  fclose(files->grid_codes_neutrophil_null_active);
  fclose(files->grid_codes_neutrophil_null_apoptotic);

  /* monocyte grid codes */
  fclose(files->grid_codes_monocyte_null_resting);
  fclose(files->grid_codes_monocyte_null_active);
  fclose(files->grid_codes_monocyte_null_apoptotic);


  /*
   * ---------------------------------------------------------------------------
   * cell counts
   * ---------------------------------------------------------------------------
   */

  /* epithelial cell counts */
  fclose(files->cell_counts_epithelial_healthy);
  fclose(files->cell_counts_epithelial_infected);
  fclose(files->cell_counts_epithelial_apoptotic);
  fclose(files->cell_counts_epithelial_phagocytosed);
  fclose(files->cell_counts_epithelial_burst);
  fclose(files->cell_counts_epithelial_ratio);

  /* infectious virus cell counts */
  fclose(files->cell_counts_virus_inf_extracellular);
  fclose(files->cell_counts_virus_inf_intracellular);
  fclose(files->cell_counts_virus_inf_total);

  /* non-infectious virus cell counts */
  fclose(files->cell_counts_virus_ninf_extracellular);
  fclose(files->cell_counts_virus_ninf_total);

  /* virus cell counts */
  fclose(files->cell_counts_virus_extracellular_total);
  fclose(files->cell_counts_virus_total);

  /* macrophage cell counts */
  fclose(files->cell_counts_macrophages_alveolar_resting);
  fclose(files->cell_counts_macrophages_alveolar_active);
  fclose(files->cell_counts_macrophages_alveolar_apoptotic);
  fclose(files->cell_counts_macrophages_alveolar);

  /* nk cell counts */
  fclose(files->cell_counts_nks_cytotoxic_resting);
  fclose(files->cell_counts_nks_cytotoxic_active);
  fclose(files->cell_counts_nks_cytotoxic_apoptotic);
  fclose(files->cell_counts_nks_cytotoxic);

  /* neutrophil cell counts */
  fclose(files->cell_counts_neutrophils_null_resting);
  fclose(files->cell_counts_neutrophils_null_active);
  fclose(files->cell_counts_neutrophils_null_apoptotic);
  fclose(files->cell_counts_neutrophils_null);

  /* monocytes cell counts */
  fclose(files->cell_counts_monocytes_null_resting);
  fclose(files->cell_counts_monocytes_null_active);
  fclose(files->cell_counts_monocytes_null_apoptotic);
  fclose(files->cell_counts_monocytes_null);


  /*
   * ---------------------------------------------------------------------------
   * Reaction diffusion systems
   * ---------------------------------------------------------------------------
   */

  /* type I interferon */
  fclose(files->rd_ifn_1);
  fclose(files->rd_total_ifn_1);

  /* il-6 */
  fclose(files->rd_il_6);
  fclose(files->rd_total_il_6);

  /* il-10 */
  fclose(files->rd_il_10);
  fclose(files->rd_total_il_10);

  /* il-12 */
  fclose(files->rd_il_12);
  fclose(files->rd_total_il_12);

  /* g_csf */
  fclose(files->rd_g_csf);
  fclose(files->rd_total_g_csf);

  /* m_csf */
  fclose(files->rd_m_csf);
  fclose(files->rd_total_m_csf);

  /* epi_inf_ck */
  fclose(files->rd_epi_inf_ck);
  fclose(files->rd_total_epi_inf_ck);

  /* NK apoptosis chemokine */
  fclose(files->rd_apop_ck);
  fclose(files->rd_total_apop_ck);

  /* macrophage phagocytosis chemokine */
  fclose(files->rd_phago_ck);
  fclose(files->rd_total_phago_ck);

  /* infectious virus */
  fclose(files->rd_virus_inf);
  fclose(files->rd_total_virus_inf);

  /* non-infectious virus */
  fclose(files->rd_virus_ninf);
  fclose(files->rd_total_virus_ninf);

  /* virus */
  fclose(files->rd_total_virus);
}
