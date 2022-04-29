#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>

#include "../cam_viralrepl.h"
#include "vr_pcv4_discon.h"
#include "../../cam_error.h"
#include "../../cam_random.h"
#include "../../cam_util.h"


/* struct containing type specific parameters */
static vr_pcv4_discon_params common_params;


/*
 * =============================================================================
 * Set the vr_pcv4_discon specific parameters
 * =============================================================================
 */
static int
vr_pcv4_discon_set_params (YAML::Node* cfg)
{
  try {
    /* dimension */
    common_params.dim = (*cfg)["viralrepl_pcv4_discon_dim"].as<int>();

    /* RK params */
    common_params.rk_step = (*cfg)["viralrepl_pcv4_discon_rk_steps"].as<int>();
    common_params.rk_step_d = (*cfg)["viralrepl_pcv4_discon_rk_steps"].as<double>();

    /* rates */
    common_params.r_U = (*cfg)["viralrepl_pcv4_discon_uncoat_rate"].as<double>();
    common_params.r_P = (*cfg)["viralrepl_pcv4_discon_rna_rate"].as<double>();
    common_params.lam_R = (*cfg)["viralrepl_pcv4_discon_rna_decay_rate"].as<double>();
    common_params.r_S = (*cfg)["viralrepl_pcv4_discon_protein_rate"].as<double>();
    common_params.r_A = (*cfg)["viralrepl_pcv4_discon_assembly_rate"].as<double>();
    common_params.lam_P = (*cfg)["viralrepl_pcv4_discon_protein_decay_rate"].as<double>();
    common_params.r_E = (*cfg)["viralrepl_pcv4_discon_export_rate"].as<double>();
    common_params.r_scale = (*cfg)["viralrepl_pcv4_discon_r_scale"].as<double>();
    common_params.A_threshold = (*cfg)["viralrepl_pcv4_discon_a_threshold"].as<double>();
    common_params.r_max = (*cfg)["viralrepl_pcv4_discon_rep_max"].as<double>();
    common_params.r_half = (*cfg)["viralrepl_pcv4_discon_rep_half"].as<double>();

    /* ifn */
    common_params.ifn_ec50 = (*cfg)["viralrepl_pcv4_discon_ifn_ec50"].as<double>();
    common_params.ifn_hill = (*cfg)["viralrepl_pcv4_discon_ifn_hill"].as<double>();

    common_params.tanh_check = true;
    common_params.A_store = 0.0;

    common_params.do_not_use_tanh_export = false;

  } catch (std::exception &e) {
    CAM_ERROR(e.what(), CAM_ERROR_IO);
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * RHS of ode system
 * =============================================================================
 */
static int
func (double t, const double y[], double f[], void* params)
{
  vr_pcv4_discon_params* par = (vr_pcv4_discon_params*) params;

  /* source */
  const double src = par->source;

  /* rates */
  const double hill = cam_util_hill_function(par->ifn_ec50, par->ifn, par->ifn_hill);

  const double rU = par->r_U;
  const double lamR = par->lam_R;
  const double rP = par->r_P * hill; //( ec50 / (ec50 + pifn) );
  const double rS = par->r_S * hill; //( ec50 / (ec50 + pifn) );
  const double rA = par->r_A;
  const double lamP = par->lam_P;
  const double rE = par->r_E;
  const double rSC = par->r_scale;
  const double aTH = par->A_threshold;
  const double repmax = par->r_max;
  const double rephalf = par->r_half;

  /* RHS: V = y[0], U = y[1], R = y[2], P = y[3], A = y[4] */
  f[0] = src - rU*y[0];
  f[1] = rU*y[0] - rP*y[1];
  f[2] = rP*y[1] + repmax * (rephalf*y[2] / (rephalf + y[2])) - lamR*y[2];
  f[3] = rS*y[2] - rA*y[3] - lamP*y[3];

  double tnh = 0.0;
  if (par->do_not_use_tanh_export) {
    tnh = 1.0;
  } else if (y[4] >= aTH) {
    par->do_not_use_tanh_export = true;
    tnh = 1.0;
  }
  f[4] = rA*y[3] - rE * tnh * y[4];

  return GSL_SUCCESS;
}


/*
 * =============================================================================
 * Jacobian for ODE system
 * =============================================================================
 */
static int
jac (double t, const double y[], double* dfdy, double dfdt[], void* params)
{
  (void) (t); /* avoid unused parameter warning */
  vr_pcv4_discon_params* par = (vr_pcv4_discon_params*) params;

  /* rates */
  const double hill = cam_util_hill_function(par->ifn_ec50, par->ifn, par->ifn_hill);

  const double rU = par->r_U;
  const double lamR = par->lam_R;
  const double rP = par->r_P * hill; //( ec50 / (ec50 + pifn) );
  const double rS = par->r_S * hill; //( ec50 / (ec50 + pifn) );
  const double rA = par->r_A;
  const double lamP = par->lam_P;
  const double rE = par->r_E;
  const double rSC = par->r_scale;
  const double aTH = par->A_threshold;
  const double repmax = par->r_max;
  const double rephalf = par->r_half;

  /* dimension */
  const int dim = par->dim;

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, dim, dim);
  gsl_matrix* m = &dfdy_mat.matrix;

  /* f[0] */
  gsl_matrix_set (m, 0, 0, -rU);
  gsl_matrix_set (m, 0, 1, 0.0);
  gsl_matrix_set (m, 0, 2, 0.0);
  gsl_matrix_set (m, 0, 3, 0.0);
  gsl_matrix_set (m, 0, 4, 0.0);

  /* f[1] */
  gsl_matrix_set (m, 1, 0, rU);
  gsl_matrix_set (m, 1, 1, -rP);
  gsl_matrix_set (m, 1, 2, 0.0);
  gsl_matrix_set (m, 1, 3, 0.0);
  gsl_matrix_set (m, 1, 4, 0.0);

  /* f[2] */
  gsl_matrix_set (m, 2, 0, 0.0);
  gsl_matrix_set (m, 2, 1, rP);
  const double dv = ((repmax*rephalf*rephalf) / ((rephalf + y[2])*(rephalf + y[2]))) - lamR;
  gsl_matrix_set (m, 2, 2, dv);
  gsl_matrix_set (m, 2, 3, 0.0);
  gsl_matrix_set (m, 2, 4, 0.0);

  /* f[3] */
  gsl_matrix_set (m, 3, 0, 0.0);
  gsl_matrix_set (m, 3, 1, 0.0);
  gsl_matrix_set (m, 3, 2, rS);
  gsl_matrix_set (m, 3, 3, -rA - lamP);
  gsl_matrix_set (m, 3, 4, 0.0);

  /* f[4] */
  double tnh = 0.0;
  if (par->do_not_use_tanh_export) {
    tnh = 1.0;
  } else if (y[4] >= aTH) {
    par->do_not_use_tanh_export = true;
    tnh = 1.0;
  }

  double tmp = -rE * tnh;
  if (y[4] == aTH) {
    tmp = GSL_NEGINF;
    return GSL_FAILURE;
  }

  gsl_matrix_set (m, 4, 0, 0.0);
  gsl_matrix_set (m, 4, 1, 0.0);
  gsl_matrix_set (m, 4, 2, 0.0);
  gsl_matrix_set (m, 4, 3, rA);
  gsl_matrix_set (m, 4, 4, tmp);

  for (int i = 0; i < dim; i++)
    dfdt[i] = 0.0;

  return GSL_SUCCESS;
}


/*
 * =============================================================================
 * Allocate the vr_pcv4_discon specific params (stored in ViralRepl object)
 * =============================================================================
 */
static void*
vr_pcv4_discon_alloc_params ()
{
  vr_pcv4_discon_params* params = (vr_pcv4_discon_params*) malloc(sizeof(vr_pcv4_discon_params));
  if (params == nullptr) {
    CAM_ERROR("cannot create vr_pcv4_discon_params", nullptr);
  }

  return params;
}


/*
 * =============================================================================
 * Allocate the vr_pcv4_discon specific data (stored in ViralRepl object)
 * =============================================================================
 */
static void*
vr_pcv4_discon_alloc_data ()
{
  vr_pcv4_discon_data* data = (vr_pcv4_discon_data*) malloc(sizeof(vr_pcv4_discon_data));
  if (data == nullptr) {
    CAM_ERROR("cannot create vr_pcv4_discon_data", nullptr);
  }

  /* GSL system */
  data->sys = (gsl_odeiv2_system*) malloc(sizeof(gsl_odeiv2_system));
  if (data == nullptr) {
    CAM_ERROR("cannot create gsl_odeiv2_system", nullptr);
  }

  return data;
}


/*
 * =============================================================================
 * Initialise the vr_pcv4_discon specific data
 * =============================================================================
 */
static int
vr_pcv4_discon_init (void* vpar, void* vdat, size_t dim)
{
  vr_pcv4_discon_params* par = (vr_pcv4_discon_params*) vpar;
  vr_pcv4_discon_data* dat = (vr_pcv4_discon_data*) vdat;

  if (common_params.dim != dim) {
    CAM_ERROR("ODE dimensions do not match", CAM_ERROR_VALUE);
  }

  /* set the local params */
  /* dimension */
  par->dim = common_params.dim;

  /* RK params */
  par->rk_step = common_params.rk_step;
  par->rk_step_d = common_params.rk_step_d;

  /* rates */
  par->r_U = common_params.r_U;
  par->r_P = common_params.r_P;
  par->lam_R = common_params.lam_R;
  par->r_S = common_params.r_S;
  par->r_A = common_params.r_A;
  par->lam_P = common_params.lam_P;
  par->r_E = common_params.r_E;
  par->r_scale = common_params.r_scale;
  par->A_threshold = common_params.A_threshold;
  par->r_max = common_params.r_max;
  par->r_half = common_params.r_half;

  /* ifn */
  par->ifn_ec50 = common_params.ifn_ec50;
  par->ifn_hill = common_params.ifn_hill;

  par->tanh_check = common_params.tanh_check;
  par->A_store = common_params.A_store;

  par->do_not_use_tanh_export = common_params.do_not_use_tanh_export;


  /* set the local data */

  /* initialise the GSL system */
  dat->sys->function = &func;
  dat->sys->jacobian = &jac;
  dat->sys->dimension = dim;
  dat->sys->params = par;

  /* allocate the GSL driver */
  dat->d = gsl_odeiv2_driver_alloc_y_new (dat->sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
  if (dat->d == nullptr) {
    CAM_ERROR("unable to create gsl_odeiv2_driver", CAM_ERROR_NOMEM);
  }

  /* initialise the solution */
  dat->sol = (double*) calloc(dim, sizeof(double));
  if (dat->sol == nullptr) {
    CAM_ERROR("unable to create gsl_odeiv2_driver", CAM_ERROR_NOMEM);
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Set the input for the ODE system (done each time step)
 * =============================================================================
 */
static void
vr_pcv4_discon_set_input (void* vpar, void* vinput, void* vifn)
{
  vr_pcv4_discon_params* par = (vr_pcv4_discon_params*) vpar;

  /* set the input data required for this ODE system */
  double* input = (double*) vinput;
  par->source = (*input);

  double* ifn = (double*) vifn;
  par->ifn = (*ifn);
}


/*
 * =============================================================================
 * Solve the ODE system
 * =============================================================================
 */
static int
vr_pcv4_discon_solve (void* vpar, void* vdat, double t_start, double t_end, void* vexp)
{
  vr_pcv4_discon_params* par = (vr_pcv4_discon_params*) vpar;
  vr_pcv4_discon_data* dat = (vr_pcv4_discon_data*) vdat;
  double* exp = (double*) vexp;

  /* check whether tanh should be used */
  double tnh = 0.0;
  if (par->do_not_use_tanh_export) {
    tnh = 1.0;
  } else if (dat->sol[4] >= par->A_threshold) {
    par->do_not_use_tanh_export = true;
    tnh = 1.0;
  }

  /* set the export */
  (*exp) = par->r_E * tnh * dat->sol[4];

  /* Set time scales */
  double t = 0.0, t1 = t_end - t_start;

  for (int i = 1; i <= par->rk_step; i++) {
    double ti = i * t1 / par->rk_step_d;

    int status = gsl_odeiv2_driver_apply (dat->d, &t, ti, dat->sol);
    if (status != GSL_SUCCESS) {
      CAM_ERROR("Error in GSL solve", CAM_ERROR_GSL);
    }
  }

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Free the vr_pcv4_discon specific data
 * =============================================================================
 */
static void
vr_pcv4_discon_free (void* vdat)
{
  vr_pcv4_discon_data* dat = (vr_pcv4_discon_data*) vdat;

  /* solution */
  free(dat->sol);

  /* GSL driver */
  gsl_odeiv2_driver_free (dat->d);

  /* GSL system */
  free(dat->sys);
}


/*
 * =============================================================================
 * Write the vr_pcv4_discon output
 * =============================================================================
 */
static void
vr_pcv4_discon_output (void* vpar, void* vdat, const double time, FILE* outf)
{
  vr_pcv4_discon_params* par = (vr_pcv4_discon_params*) vpar;
  vr_pcv4_discon_data* dat = (vr_pcv4_discon_data*) vdat;

  double tmp_tanh = 1.0;
  if (par->tanh_check) {
    tmp_tanh = 0.5 * (1.0 + tanh(par->r_scale*(dat->sol[4] - par->A_threshold)));
  }
  double exp = par->r_E * tmp_tanh * dat->sol[4];

  fprintf(outf, "%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", time, dat->sol[0], dat->sol[1], dat->sol[2], dat->sol[3], dat->sol[4], exp);
}


/*
 * =============================================================================
 * Get the solution data at specified index
 * =============================================================================
 */
static double
vr_pcv4_discon_get_sol (void* vdat, const int index)
{
  vr_pcv4_discon_data* dat = (vr_pcv4_discon_data*) vdat;

  return (dat->sol[index]);
}


/*
 * =============================================================================
 * Get whether we are exporting or not
 * =============================================================================
 */
static bool
vr_pcv4_discon_get_exporting (void* vpar)
{
  vr_pcv4_discon_params* par = (vr_pcv4_discon_params*) vpar;

  return (par->do_not_use_tanh_export);
}


/*
 * =============================================================================
 * Statically define the pcv4_ode type
 * =============================================================================
 */
static const cam_vr_type vr_pcv4_discon_type =
{
  "pcv4_discon",
  5,
  &vr_pcv4_discon_set_params,
  &vr_pcv4_discon_alloc_params,
  &vr_pcv4_discon_alloc_data,
  &vr_pcv4_discon_init,
  &vr_pcv4_discon_set_input,
  &vr_pcv4_discon_solve,
  &vr_pcv4_discon_free,
  &vr_pcv4_discon_output,
  &vr_pcv4_discon_get_sol,
  &vr_pcv4_discon_get_exporting
};


/*
 * =============================================================================
 * Set the global vr_pcv4_discon type
 * =============================================================================
 */
const cam_vr_type* cam_vr_type_pcv4_discon = &vr_pcv4_discon_type;
