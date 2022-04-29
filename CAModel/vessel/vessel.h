#ifndef __VESSEL_H__
#define __VESSEL_H__

/* vessel specific parameters */
typedef struct
{
  /* vessel radius */
  double vessel_radius;
  
  /* cytokine uptake by vessels */
  double vessel_blood_uptake_ifn_1;
  double vessel_lymph_uptake_ifn_1;

  double vessel_blood_uptake_vir_inf;
  double vessel_lymph_uptake_vir_inf;

} vessel_params;

#endif /* __EPITHELIAL_H__ */
