#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

#include "../cam_recepexpr.h"
#include "re_pcv2_discon.h"

#include "../../cam_error.h"
#include "../../cam_util.h"


/* struct containing type specific parameters */
static re_pcv2_discon_params common_params;


/*
 * =============================================================================
 * Set the re_pcv2_discon specific parameters
 * =============================================================================
 */
static int
re_pcv2_discon_set_params (YAML::Node* cfg)
{
  try {
    /* dimension */
    common_params.dim = (*cfg)["recepexpr_pcv2_discon_dim"].as<int>();

    /* RK params */
    common_params.rk_step = (*cfg)["recepexpr_pcv2_discon_rk_steps"].as<int>();
    common_params.rk_step_d = (*cfg)["recepexpr_pcv2_discon_rk_steps"].as<double>();

    /* number of receptors */
    common_params.r_tot_min = (*cfg)["recepexpr_pcv2_discon_r_tot_min"].as<int>();
    common_params.r_tot_max = (*cfg)["recepexpr_pcv2_discon_r_tot_max"].as<int>();
    common_params.r_tot = 0;

    /* rates */
    common_params.r_bind = (*cfg)["recepexpr_pcv2_discon_bind_rate"].as<double>();
    common_params.r_recycle = (*cfg)["recepexpr_pcv2_discon_recycle_rate"].as<double>();
    common_params.r_endo = (*cfg)["recepexpr_pcv2_discon_endocytosis_rate"].as<double>();
    common_params.r_release = (*cfg)["recepexpr_pcv2_discon_release_rate"].as<double>();

    common_params.r_scale = (*cfg)["recepexpr_pcv2_discon_r_scale"].as<double>();
    common_params.nvir_threshold = (*cfg)["recepexpr_pcv2_discon_nvir_threshold"].as<double>();

    /* ifn */
    common_params.ifn_ec50 = (*cfg)["recepexpr_pcv2_discon_ifn_ec50"].as<double>();
    common_params.ifn_hill = (*cfg)["recepexpr_pcv2_discon_ifn_hill"].as<double>();

    common_params.tanh_unsafe = false;
    common_params.nvir_store = 0.0;

    common_params.do_not_use_tanh_binding = false;

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
  (void) (t); /* avoid unused parameter warning */
  re_pcv2_discon_params* par = (re_pcv2_discon_params*) params;

  /* input */
  const double n_v = par->num_virs;

  /* rates */
  const double hill = cam_util_hill_function(par->ifn_ec50, par->ifn, par->ifn_hill);

  double tnh = 0.0;
  if (par->do_not_use_tanh_binding) {
    tnh = 1.0;
  } else if (n_v >= par->nvir_threshold) {
    par->do_not_use_tanh_binding = true;
    tnh = 1.0;
  }

  const double r_b = par->r_bind * hill * tnh;
  const double r_rc = par->r_recycle;
  const double r_e = par->r_endo * hill; //( ec50 / (ec50 + pifn) );
  const double r_re = par->r_release;

  /* RHS: Reu = y[0], Reb = y[1], Rib = y[2], Riu = y[3] */
  f[0] = -r_b*n_v*y[0] + r_rc*y[3];
  f[1] = r_b*n_v*y[0] - r_e*y[1];
  f[2] = r_e*y[1] - r_re*y[2];
  f[3] = r_re*y[2] - r_rc*y[3];

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
  (void) (y); /* avoid unused parameter warning */
  re_pcv2_discon_params* par = (re_pcv2_discon_params*) params;

  /* input */
  const double n_v = par->num_virs;

  /* rates */
  const double hill = cam_util_hill_function(par->ifn_ec50, par->ifn, par->ifn_hill);

  double tnh = 0.0;
  if (par->do_not_use_tanh_binding) {
    tnh = 1.0;
  } else if (n_v >= par->nvir_threshold) {
    par->do_not_use_tanh_binding = true;
    tnh = 1.0;
  }

  const double r_b = par->r_bind * hill * tnh;
  const double r_rc = par->r_recycle;
  const double r_e = par->r_endo * hill; //( ec50 / (ec50 + pifn) );
  const double r_re = par->r_release;

  /* dimension */
  const int dim = par->dim;

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, dim, dim);
  gsl_matrix* m = &dfdy_mat.matrix;

  /* f[0] = -r_b*n_v*y[0] + r_rc*y[3]; */
  gsl_matrix_set (m, 0, 0, -r_b*n_v);
  gsl_matrix_set (m, 0, 1, 0.0);
  gsl_matrix_set (m, 0, 2, 0.0);
  gsl_matrix_set (m, 0, 3, r_rc);

  /* f[1] = r_b*n_v*y[0] - r_e*y[1]; */
  gsl_matrix_set (m, 1, 0, r_b*n_v);
  gsl_matrix_set (m, 1, 1, -r_e);
  gsl_matrix_set (m, 1, 2, 0.0);
  gsl_matrix_set (m, 1, 3, 0.0);

  /* f[2] = r_e*y[1] - r_re*y[2]; */
  gsl_matrix_set (m, 2, 0, 0.0);
  gsl_matrix_set (m, 2, 1, r_e);
  gsl_matrix_set (m, 2, 2, -r_re);
  gsl_matrix_set (m, 2, 3, 0.0);

  /* f[3] = r_re*y[2] - r_rc*y[3]; */
  gsl_matrix_set (m, 3, 0, 0.0);
  gsl_matrix_set (m, 3, 1, 0.0);
  gsl_matrix_set (m, 3, 2, r_re);
  gsl_matrix_set (m, 3, 3, -r_rc);

  for (int i = 0; i < dim; i++)
    dfdt[i] = 0.0;

  return GSL_SUCCESS;
}


/*
 * =============================================================================
 * Allocate the re_pcv2_discon specific params (stored in RecepExpr object)
 * =============================================================================
 */
static void*
re_pcv2_discon_alloc_params ()
{
  re_pcv2_discon_params* params = (re_pcv2_discon_params*) malloc(sizeof(re_pcv2_discon_params));
  if (params == nullptr) {
    CAM_ERROR("cannot create re_pcv2_discon_params", nullptr);
  }

  return params;
}


/*
 * =============================================================================
 * Allocate the re_pcv2_discon specific data (stored in RecepExpr object)
 * =============================================================================
 */
static void*
re_pcv2_discon_alloc_data ()
{
  re_pcv2_discon_data* data = (re_pcv2_discon_data*) malloc(sizeof(re_pcv2_discon_data));
  if (data == nullptr) {
    CAM_ERROR("cannot create re_pcv2_discon_data", nullptr);
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
 * Initialise the re_pcv2_discon specific data
 * =============================================================================
 */
static int
re_pcv2_discon_init (void* vpar, void* vdat, size_t dim)
{
  re_pcv2_discon_params* par = (re_pcv2_discon_params*) vpar;
  re_pcv2_discon_data* dat = (re_pcv2_discon_data*) vdat;

  if (common_params.dim != dim) {
    CAM_ERROR("ODE dimensions do not match", CAM_ERROR_VALUE);
  }

  /* initialise the local parameters */
  /* dimension */
  par->dim = common_params.dim;

  /* RK params */
  par->rk_step = common_params.rk_step;
  par->rk_step_d = common_params.rk_step_d;

  /* number of receptors */
  par->r_tot_min = common_params.r_tot_min;
  par->r_tot_max = common_params.r_tot_max;
  par->r_tot = cam_util_choose_random_index(common_params.r_tot_min, common_params.r_tot_max);

  /* rates */
  par->r_bind = common_params.r_bind;
  par->r_recycle = common_params.r_recycle;
  par->r_endo = common_params.r_endo;
  par->r_release = common_params.r_release;

  par->r_scale = common_params.r_scale;
  par->nvir_threshold = common_params.nvir_threshold;

  /* ifn */
  par->ifn_ec50 = common_params.ifn_ec50;
  par->ifn_hill = common_params.ifn_hill;

  par->tanh_unsafe = common_params.tanh_unsafe;
  par->nvir_store = common_params.nvir_store;

  par->do_not_use_tanh_binding = common_params.do_not_use_tanh_binding;


  /* initialise the local data */

  /* initialise the GSL system */
  dat->sys->function = &func;
  dat->sys->jacobian = &jac;
  dat->sys->dimension = dim;
  dat->sys->params = par;

  /* allocate the GSL driver */
  dat->d = gsl_odeiv2_driver_alloc_y_new(dat->sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
  if (dat->d == nullptr) {
    CAM_ERROR("unable to create gsl_odeiv2_driver", CAM_ERROR_NOMEM);
  }

  /* initialise the solution */
  dat->sol = (double*) calloc(dim, sizeof(double));
  if (dat->sol == nullptr) {
    CAM_ERROR("unable to create gsl_odeiv2_driver", CAM_ERROR_NOMEM);
  }

  /* set the initial conditions */
  dat->sol[0] = par->r_tot;

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * Set the input for the ODE system (done each time step)
 * =============================================================================
 */
static void
re_pcv2_discon_set_input (void* vpar, void* vinput, void* vifn)
{
  re_pcv2_discon_params* par = (re_pcv2_discon_params*) vpar;

  /* set the input data required for this ODE system */
  double* input = (double*) vinput;
  par->num_virs = (*input);

  double* ifn = (double*) vifn;
  par->ifn = (*ifn);
}


/*
 * =============================================================================
 * Solve the ODE system
 * =============================================================================
 */
static int
re_pcv2_discon_solve (void* vpar, void* vdat, double t_start, double t_end, void* vupt, void* vsrc)
{
  re_pcv2_discon_params* par = (re_pcv2_discon_params*) vpar;
  re_pcv2_discon_data* dat = (re_pcv2_discon_data*) vdat;
  double* uptake = (double*) vupt;
  double* source = (double*) vsrc;

  /* set the uptake */
  const double hill = cam_util_hill_function(par->ifn_ec50, par->ifn, par->ifn_hill);

  double tnh = 0.0;
  if (par->do_not_use_tanh_binding) {
    tnh = 1.0;
  } else if (par->num_virs >= par->nvir_threshold) {
    par->do_not_use_tanh_binding = true;
    tnh = 1.0;
  }
  (*uptake) = par->r_bind * hill * tnh * dat->sol[0];

  /* set the source */
  (*source) = par->r_release * dat->sol[2]; // r_release * R_ib

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
 * Free the re_pcv2_discon specific data
 * =============================================================================
 */
static void
re_pcv2_discon_free (void* vdat)
{
  re_pcv2_discon_data* dat = (re_pcv2_discon_data*) vdat;

  /* solution */
  free(dat->sol);

  /* GSL driver */
  gsl_odeiv2_driver_free (dat->d);

  /* GSL system */
  free(dat->sys);
}


/*
 * =============================================================================
 * Write the re_pcv2_discon output
 * =============================================================================
 */
static void
re_pcv2_discon_output (void* vpar, void* vdat, const double time, FILE* outf)
{
  re_pcv2_discon_params* par = (re_pcv2_discon_params*) vpar;
  re_pcv2_discon_data* dat = (re_pcv2_discon_data*) vdat;
  double exp = par->r_release * dat->sol[2]; // r_release * R_ib

  fprintf(outf, "%.5e %.5e %.5e %.5e %.5e %.5e\n", time, dat->sol[0], dat->sol[1], dat->sol[2], dat->sol[3], exp);
}


/*
 * =============================================================================
 * Get the solution data at specified index
 * =============================================================================
 */
static double
re_pcv2_discon_get_sol (void* vdat, const int index)
{
  re_pcv2_discon_data* dat = (re_pcv2_discon_data*) vdat;

  return (dat->sol[index]);
}


/*
 * =============================================================================
 * Statically define the re_pcv2_discon type
 * =============================================================================
 */
static const cam_re_type re_pcv2_discon_type =
{
  "pcv2_discon",
  4,
  &re_pcv2_discon_set_params,
  &re_pcv2_discon_alloc_params,
  &re_pcv2_discon_alloc_data,
  &re_pcv2_discon_init,
  &re_pcv2_discon_set_input,
  &re_pcv2_discon_solve,
  &re_pcv2_discon_free,
  &re_pcv2_discon_output,
  &re_pcv2_discon_get_sol
};


/*
 * =============================================================================
 * Set the global re_pcv2_discon type
 * =============================================================================
 */
const cam_re_type* cam_re_type_pcv2_discon = &re_pcv2_discon_type;
