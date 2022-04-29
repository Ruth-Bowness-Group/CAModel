#ifndef __CAM_OUTPUT_H__
#define __CAM_OUTPUT_H__

#include <cstdio>
#include <string>


/*******************************************************************************
 * Simple macros for simplifying output
 ******************************************************************************/

/* reset the grid codes */
#define CAM_OUTPUT_RESET_GRID_CODES(val) \
  do { \
    for (int x = 0; x < grid_size; x++) { \
      for (int y = 0; y < grid_size; y++) { \
        grid_codes[x][y] = (val); \
      } \
    } \
  } while (0)

/* write the grid codes */
#define CAM_OUTPUT_WRITE_GRID_CODES(fp) \
  do { \
    for (int x = 0; x < grid_size; x++) { \
      for (int y = 0; y < grid_size; y++) { \
        fprintf((fp), "%d\t", grid_codes[x][y]); \
      } \
    } \
    fprintf((fp), "\n"); \
  } while (0)


/*******************************************************************************
 * Output file structure
 ******************************************************************************/

struct cam_output_files {
  /* constants */
  FILE* consts;

  /*
   * ---------------------------------------------------------------------------
   * grid output codes
   * ---------------------------------------------------------------------------
   */

  /* location grid codes */
  FILE* grid_codes;

  /* macrophage grid codes */
  FILE* grid_codes_macrophage_alveolar_resting;
  FILE* grid_codes_macrophage_alveolar_active;
  FILE* grid_codes_macrophage_alveolar_apoptotic;

  /* nk grid codes */
  FILE* grid_codes_nk_cytotoxic_resting;
  FILE* grid_codes_nk_cytotoxic_active;
  FILE* grid_codes_nk_cytotoxic_apoptotic;

  /* neutrophil grid codes */
  FILE* grid_codes_neutrophil_null_resting;
  FILE* grid_codes_neutrophil_null_active;
  FILE* grid_codes_neutrophil_null_apoptotic;

  /* monocyte grid codes */
  FILE* grid_codes_monocyte_null_resting;
  FILE* grid_codes_monocyte_null_active;
  FILE* grid_codes_monocyte_null_apoptotic;


  /*
   * ---------------------------------------------------------------------------
   * cell counts
   * ---------------------------------------------------------------------------
   */

  /* epithelial cell counts */
  FILE* cell_counts_epithelial_healthy;
  FILE* cell_counts_epithelial_infected;
  FILE* cell_counts_epithelial_apoptotic;
  FILE* cell_counts_epithelial_phagocytosed;
  FILE* cell_counts_epithelial_burst;
  FILE* cell_counts_epithelial_ratio;

  /* infectious virus cell counts */
  FILE* cell_counts_virus_inf_extracellular;
  FILE* cell_counts_virus_inf_intracellular;
  FILE* cell_counts_virus_inf_total;

  /* non-infectious virus cell counts */
  FILE* cell_counts_virus_ninf_extracellular;
  FILE* cell_counts_virus_ninf_total;

  /* virus cell counts */
  FILE* cell_counts_virus_extracellular_total;
  FILE* cell_counts_virus_total;

  /* macrophage cell counts */
  FILE* cell_counts_macrophages_alveolar_resting;
  FILE* cell_counts_macrophages_alveolar_active;
  FILE* cell_counts_macrophages_alveolar_apoptotic;
  FILE* cell_counts_macrophages_alveolar;

  /* nk cell counts */
  FILE* cell_counts_nks_cytotoxic_resting;
  FILE* cell_counts_nks_cytotoxic_active;
  FILE* cell_counts_nks_cytotoxic_apoptotic;
  FILE* cell_counts_nks_cytotoxic;

  /* neutrophil cell counts */
  FILE* cell_counts_neutrophils_null_resting;
  FILE* cell_counts_neutrophils_null_active;
  FILE* cell_counts_neutrophils_null_apoptotic;
  FILE* cell_counts_neutrophils_null;

  /* monocyte cell counts */
  FILE* cell_counts_monocytes_null_resting;
  FILE* cell_counts_monocytes_null_active;
  FILE* cell_counts_monocytes_null_apoptotic;
  FILE* cell_counts_monocytes_null;


  /*
   * ---------------------------------------------------------------------------
   * reaction diffusion systems
   * ---------------------------------------------------------------------------
   */

  /* interferon */
  FILE* rd_ifn_1;
  FILE* rd_total_ifn_1;

  /* il */
  FILE* rd_il_6;
  FILE* rd_total_il_6;
  FILE* rd_il_10;
  FILE* rd_total_il_10;
  FILE* rd_il_12;
  FILE* rd_total_il_12;

  /* g_csf */
  FILE* rd_g_csf;
  FILE* rd_total_g_csf;

  /* m_csf */
  FILE* rd_m_csf;
  FILE* rd_total_m_csf;

  /* epi_inf_ck */
  FILE* rd_epi_inf_ck;
  FILE* rd_total_epi_inf_ck;

  /* NK apoptosis chemokine */
  FILE* rd_apop_ck;
  FILE* rd_total_apop_ck;

  /* macrophage phagocytosis chemokine */
  FILE* rd_phago_ck;
  FILE* rd_total_phago_ck;

  /* infectious virus */
  FILE* rd_virus_inf;
  FILE* rd_total_virus_inf;

  /* non-infectious virus */
  FILE* rd_virus_ninf;
  FILE* rd_total_virus_ninf;

  /* virus */
  FILE* rd_total_virus;
};

/* open the output files */
int cam_output_open(const std::string output_folder, struct cam_output_files* files);

/* close the output files */
void cam_output_close(struct cam_output_files* files);

#endif /* __CAM_OUTPUT_H__ */
