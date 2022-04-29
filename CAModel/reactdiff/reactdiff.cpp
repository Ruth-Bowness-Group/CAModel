#include <cmath>

#include "cam_reactdiff.h"
#include "reactdiff.h"
#include "../cam_error.h"


/* struct containing reaction diffusion specific parameters */
static struct rd_params params;


/*
 * =============================================================================
 * Set the reaction diffusion specific parameters
 * =============================================================================
 */
int
cam_rd_set_params (YAML::Node* cfg)
{
  try {
    /* Spatial step of the grid */
    params.rd_spatial_step = (*cfg)["rd_spatial_step"].as<double>();

    /* T1 IFN */
    params.rd_ifn_1_diffusion = (*cfg)["rd_ifn_1_diffusion"].as<double>();
    params.rd_ifn_1_decay = (*cfg)["rd_ifn_1_decay"].as<double>();

    /* IL-6 */
    params.rd_il_6_diffusion = (*cfg)["rd_il_6_diffusion"].as<double>();
    params.rd_il_6_decay = (*cfg)["rd_il_6_decay"].as<double>();

    /* IL-10 */
    params.rd_il_10_diffusion = (*cfg)["rd_il_10_diffusion"].as<double>();
    params.rd_il_10_decay = (*cfg)["rd_il_10_decay"].as<double>();

    /* IL-12 */
    params.rd_il_12_diffusion = (*cfg)["rd_il_12_diffusion"].as<double>();
    params.rd_il_12_decay = (*cfg)["rd_il_12_decay"].as<double>();

    /* G-CSF */
    params.rd_g_csf_diffusion = (*cfg)["rd_g_csf_diffusion"].as<double>();
    params.rd_g_csf_decay = (*cfg)["rd_g_csf_decay"].as<double>();

    /* M-CSF */
    params.rd_m_csf_diffusion = (*cfg)["rd_m_csf_diffusion"].as<double>();
    params.rd_m_csf_decay = (*cfg)["rd_m_csf_decay"].as<double>();

    /* infected epithelial chemokine */
    params.rd_epi_inf_ck_diffusion = (*cfg)["rd_epi_inf_ck_diffusion"].as<double>();
    params.rd_epi_inf_ck_decay = (*cfg)["rd_epi_inf_ck_decay"].as<double>();

    /* apoptosis chemokine */
    params.rd_apop_ck_diffusion = (*cfg)["rd_apop_ck_diffusion"].as<double>();
    params.rd_apop_ck_decay = (*cfg)["rd_apop_ck_decay"].as<double>();

    /* phagocytosis chemokine */
    params.rd_phago_ck_diffusion = (*cfg)["rd_phago_ck_diffusion"].as<double>();
    params.rd_phago_ck_decay = (*cfg)["rd_phago_ck_decay"].as<double>();

    /* infectious virus */
    params.rd_virus_inf_diffusion = (*cfg)["rd_virus_inf_diffusion"].as<double>();
    params.rd_virus_inf_decay = (*cfg)["rd_virus_inf_decay"].as<double>();

    /* non-infectious virus */
    params.rd_virus_ninf_diffusion = (*cfg)["rd_virus_ninf_diffusion"].as<double>();
    params.rd_virus_ninf_decay = (*cfg)["rd_virus_ninf_decay"].as<double>();

  } catch (std::exception &e) {
    CAM_ERROR(e.what(), CAM_ERROR_IO);
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * ReactDiff constuctor
 * =============================================================================
 */
ReactDiff::ReactDiff (RDType tp, int* status)
{
  double df = 0.0;
  double dc = 0.0;
  double init_val = 0.0;

  /* Initialise the type of RD system */
  set_type(tp);

  /* Depending on type, set the diffusion constant */
  switch(tp) {
    case RDType::ifn_1:
      df = params.rd_ifn_1_diffusion;
      dc = params.rd_ifn_1_decay;
      break;
    case RDType::il_6:
      df = params.rd_il_6_diffusion;
      dc = params.rd_il_6_decay;
      break;
    case RDType::il_10:
      df = params.rd_il_10_diffusion;
      dc = params.rd_il_10_decay;
      break;
    case RDType::il_12:
      df = params.rd_il_12_diffusion;
      dc = params.rd_il_12_decay;
      break;
    case RDType::g_csf:
      df = params.rd_g_csf_diffusion;
      dc = params.rd_g_csf_decay;
      break;
    case RDType::m_csf:
      df = params.rd_m_csf_diffusion;
      dc = params.rd_m_csf_decay;
      break;
    case RDType::epi_inf_ck:
      df = params.rd_epi_inf_ck_diffusion;
      dc = params.rd_epi_inf_ck_decay;
      break;
    case RDType::apop_ck:
      df = params.rd_apop_ck_diffusion;
      dc = params.rd_apop_ck_decay;
      break;
    case RDType::phago_ck:
      df = params.rd_phago_ck_diffusion;
      dc = params.rd_phago_ck_decay;
      init_val = (validation_chemotaxis_and_recruitment) ? 1.0 : 0.0;
      break;
    case RDType::virus_inf:
      df = params.rd_virus_inf_diffusion;
      dc = params.rd_virus_inf_decay;
      break;
    case RDType::virus_ninf:
      df = params.rd_virus_ninf_diffusion;
      dc = params.rd_virus_ninf_decay;
      break;
    default:
      (*status) = CAM_ERROR_BADARGS;
      CAM_ERROR_VOID("unknown reaction diffusion type");
  }

  /* Set the initial conditions */
  for (int x = 0; x < grid_size; x++) {
    for (int y = 0; y < grid_size; y++) {
      diffusion[x][y] = df;
      decay[x][y] = dc;

      conc[x][y] = 0.0;

      basal_source[x][y] = (x == 40 && y == 80) ? init_val : 0.0; // set by cell
      source[x][y] = 0.0; // set by cell
      burst_source[x][y] = 0.0; // set by burst cell

      basal_uptake[x][y] = 0.0; // set by cell
      uptake[x][y] = 0.0; // set by cell
    }
  }
  max_conc = 0.0;

  /* successful */
  (*status) = CAM_SUCCESS;
}


/*
 * =============================================================================
 * Destructor
 * =============================================================================
 */
ReactDiff::~ReactDiff ()
{
  // nothing to do
}


/*
 * =============================================================================
 * Advance the reaction diffusion system to next time step
 * =============================================================================
 */
int
ReactDiff::advance ()
{
  int x = 0, y = 0;
  double tmp_ch[grid_size][grid_size] = {{0.0}};

  for (x = 0; x < grid_size; x++) {
    for (y = 0; y < grid_size; y++) {
      /* Get the cytokine values at the centre, north, east, south, west neighbours */
      const double val_c = conc[x][y];

      const double val_n = ( (y == grid_size-1) ? conc[x][y-1] : conc[x][y+1] );
      const double val_s = ( (y == 0) ? conc[x][y+1] : conc[x][y-1] );
      const double val_e = ( (x == grid_size-1) ? conc[x-1][y] : conc[x+1][y] );
      const double val_w = ( (x == 0) ? conc[x+1][y] : conc[x-1][y] );


      /* Get the diffusion values at the centre, north, east, south, west neighbours */
      const double df_c = diffusion[x][y];

      const double df_n = ( (y == grid_size-1) ? diffusion[x][y-1] : diffusion[x][y+1] );
      const double df_s = ( (y == 0) ? diffusion[x][y+1] : diffusion[x][y-1] );
      const double df_e = ( (x == grid_size-1) ? diffusion[x-1][y] : diffusion[x+1][y] );
      const double df_w = ( (x == 0) ? diffusion[x+1][y] : diffusion[x-1][y] );


      /* source, decay and uptake values */
      const double ch_src = basal_source[x][y] + source[x][y]*val_c + burst_source[x][y];
      const double ch_up = basal_uptake[x][y] + uptake[x][y]*val_c;
      const double ch_dc = decay[x][y]*val_c;


      /* Calculate diffusion values using determined neighbours */
      double ch_df;
      if ((x > 0) && (x < grid_size-1) && (y > 0) && (y < grid_size-1)) {
        /* Interior */
        ch_df = calc_diffusion(val_c, val_n, val_e, val_s, val_w, df_c, df_n, df_e, df_s, df_w);
      } else {
        /* Boundary */
        ch_df = calc_diffusion_bdry(val_c, val_n, val_e, val_s, val_w, df_c);
      }


      /* CONCENTRATION: Current value + timestep*(diffusion + source - uptake - decay) */
      tmp_ch[x][y] = val_c + timestep * (ch_df + ch_src - ch_up - ch_dc);
      if (tmp_ch[x][y] != tmp_ch[x][y]) {
        fprintf(stderr, "\trd: (type, %d)\n", static_cast<int>(type));
        fprintf(stderr, "\trd: (%d, %d)\n", x, y);
        fprintf(stderr, "\trd: \ttmp_ch[%d][%d] = %g\n", x, y, tmp_ch[x][y]);
        fprintf(stderr, "\trd: \tval_c = %g\n", val_c);
        fprintf(stderr, "\trd: \ttstep = %g\n", timestep);
        fprintf(stderr, "\trd: \tdiffusion = %g\n", ch_df);
        fprintf(stderr, "\trd: \tsource = %g\n", ch_src);
        fprintf(stderr, "\trd: \tuptake[%d][%d] = %g\n", x, y, uptake[x][y]);
        fprintf(stderr, "\trd: \tbasal_uptake[%d][%d] = %g\n", x, y, basal_uptake[x][y]);
        fprintf(stderr, "\trd: \tdecay = %g\n", ch_dc);
        CAM_ERROR("rd conc is NaN", CAM_ERROR_VALUE);
      }
      if (tmp_ch[x][y] < NEG_ZERO) {
        if (abs(tmp_ch[x][y]) < CAM_TOL) {
          tmp_ch[x][y] = 0.0;
        } else {
          fprintf(stderr, "\trd: (type, %d)\n", static_cast<int>(type));
          fprintf(stderr, "\trd: (%d, %d)\n", x, y);
          fprintf(stderr, "\trd: \ttmp_ch[%d][%d] = %g\n", x, y, tmp_ch[x][y]);
          fprintf(stderr, "\trd: \tval_c = %g\n", val_c);
          fprintf(stderr, "\trd: \ttstep = %g\n", timestep);
          fprintf(stderr, "\trd: \tdiffusion = %g\n", ch_df);
          fprintf(stderr, "\trd: \tsource = %g\n", ch_src);
          fprintf(stderr, "\trd: \tuptake = %g\n", ch_up);
          fprintf(stderr, "\trd: \tbasal_uptake[%d][%d] = %g\n", x, y, basal_uptake[x][y]);
          fprintf(stderr, "\trd: \tdecay = %g\n", ch_dc);
          CAM_ERROR("rd conc less than zero", CAM_ERROR_VALUE);
        }
      }
    }
  }

  /* update the concentration */
  max_conc = 0.0;
  for (x = 0; x < grid_size; x++) {
    for (y = 0; y < grid_size; y++) {
      conc[x][y] = ((tmp_ch[x][y] < CAM_TOL) ? 0.0 : tmp_ch[x][y]);

      if (conc[x][y] > max_conc) {
        max_conc = conc[x][y];
      }

      /* decrease the burst source */
      if (burst_source[x][y] > CAM_TOL) {
        dec_burst_source(x, y, timestep*burst_source[x][y]);
      }
    }
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Calculate the diffusion term
 * =============================================================================
 */

/* interior diffusion calculation */
double
ReactDiff::calc_diffusion (double v_c, double v_n, double v_e, double v_s, double v_w,
                           double d_c, double d_n, double d_e, double d_s, double d_w)
{
  double tmp = 0.0;
  double tmp2 = 0.0;
  const double h2 = params.rd_spatial_step * params.rd_spatial_step;

  tmp = ((d_c + d_e) / 2.0) * (v_e - v_c);
  tmp -= ((d_w + d_c) / 2.0) * (v_c - v_w);
  tmp /= h2;

  tmp2 = ((d_c + d_n) / 2.0) * (v_n - v_c);
  tmp2 -= ((d_s + d_c) / 2.0) * (v_c - v_s);
  tmp2 /= h2;

  return tmp + tmp2;
}


/* boundary diffusion calculation */
double
ReactDiff::calc_diffusion_bdry (double v_c, double v_n, double v_e, double v_s, double v_w, double d_c)
{
  const double h2 = params.rd_spatial_step * params.rd_spatial_step;

  const double tmp = v_w - 2*v_c + v_e;
  const double tmp2 = v_n - 2*v_c + v_s;

  return (d_c * (tmp + tmp2) / h2);
}
