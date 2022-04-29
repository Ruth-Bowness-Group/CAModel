#ifndef __VR_PCV4_DISCON_H__
#define __VR_PCV4_DISCON_H__

#include <gsl/gsl_odeiv2.h>

/* pcv4_discon specific parameters */
typedef struct {
  size_t dim;

  /* RK step */
  int rk_step;
  double rk_step_d;

  /* source */
  double source;

  /* rates */
  double r_U;
  double r_P;
  double lam_R;
  double r_S;
  double r_A;
  double lam_P;
  double r_E;
  double r_max;
  double r_half;

  /* ifn */
  double ifn;
  double ifn_ec50;
  double ifn_hill;

  double r_scale;
  double A_threshold;

  bool tanh_check;
  double A_store;

  bool do_not_use_tanh_export;

} vr_pcv4_discon_params;

/* type specific data */
typedef struct {
  /* GSL components */
  gsl_odeiv2_system* sys; // GSL system
  gsl_odeiv2_driver* d;  // GSL driver

  /* solution store */
  double* sol;

} vr_pcv4_discon_data;

#endif /* __VR_PCV4_DISCON_H__ */
