#ifndef __ENVIRONMENT_H__
#define __ENVIRONMENT_H__

/* struct containing environment specific params */
struct env_params
{
  /* Blood vessels */
  bool env_fixed_blood_vessel_placement;
  const std::vector<std::pair<int, int>> fixed_blood_vessels = {{grid_size/2, grid_size/2}};
  int env_initial_number_blood_vessels;

  /* Lymphatic vessels */
  bool env_fixed_lymph_vessel_placement;
  const std::vector<std::pair<int,int>> fixed_lymph_vessels = {{25, 50}, {50, 75}, {75, 50}, {50, 25}};
  int env_initial_number_lymph_vessels;

  /* Virions */
  bool env_fixed_virion_placement;
  const std::vector<std::pair<int, int>> fixed_virions = {
    {grid_size/2, grid_size/2}
  };
  int env_initial_number_virions;
  double env_initial_moi;

  /* Macrophages */
  bool env_fixed_macrophage_placement;
  const std::vector<std::pair<int, int>> fixed_macrophages = {{grid_size/2, grid_size/2}};
  int env_initial_number_macrophages;

  /* NKs */
  bool env_fixed_nk_placement;
  const std::vector<std::pair<int, int>> fixed_nks = {{65, 65}};
  int env_initial_number_nks;

  /* Neutrophils */
  bool env_fixed_neutrophil_placement;
  const std::vector<std::pair<int, int>> fixed_neutrophils = {{65, 65}};
  int env_initial_number_neutrophils;

  /* Monocytes */
  bool env_fixed_monocyte_placement;
  const std::vector<std::pair<int, int>> fixed_monocytes = {{65, 65}};
  int env_initial_number_monocytes;

  /* maximum neighbourhood */
  int env_neighbourhood_max_level;
};

#endif /* __ENVIRONMENT_H__ */
