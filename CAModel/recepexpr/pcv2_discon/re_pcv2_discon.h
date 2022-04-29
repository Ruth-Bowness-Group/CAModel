#ifndef __RE_PCV2_DISCON_H__
#define __RE_PCV2_DISCON_H__

#include <gsl/gsl_odeiv2.h>

/* pcv2_discon specific parameters */
typedef struct {
  size_t dim;

  /* RK step */
  int rk_step;
  double rk_step_d;

  /* number of receptors */
  int r_tot_min;
  int r_tot_max;
  int r_tot;

  /* number of virions */
  double num_virs;

  /* rates */
  double r_bind;
  double r_recycle;
  double r_endo;
  double r_release;

  /* ifn */
  double ifn;
  double ifn_ec50;
  double ifn_hill;

  double r_scale;
  double nvir_threshold;

  bool tanh_unsafe;
  double nvir_store;

  bool do_not_use_tanh_binding;

} re_pcv2_discon_params;

/* type specific data */
typedef struct {
  /* GSL components */
  gsl_odeiv2_system* sys; // GSL system
  gsl_odeiv2_driver* d;  // GSL driver

  /* solution store */
  double* sol;

} re_pcv2_discon_data;

#endif /* __RE_PCV2_DISCON_H__ */
