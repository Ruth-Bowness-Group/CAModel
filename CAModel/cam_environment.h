#ifndef __CAM_ENVIRONMENT_H__
#define __CAM_ENVIRONMENT_H__

#include <yaml-cpp/yaml.h>
#include <vector>

#include "cam_parameters.h"
#include "cam_location.h"
#include "cam_output.h"
#include "reactdiff/cam_reactdiff.h"
#include "vessel/cam_vessel.h"
#include "epithelial/cam_epithelial.h"


/* set the environment specific parameters */
int cam_env_set_params (YAML::Node* cfg);


/*
 * =============================================================================
 * Tissue environment class
 * =============================================================================
 */
class Environment
{
private:
  int number_of_outputs;

  /* grid of locations */
  Location* grid[grid_size][grid_size];

  /* struct containing eaction diffusion equations */
  struct cam_rd_sys rds;

  /* blood and lymphatic vessels */
  std::vector<Vessel*> vessels;
  int init_vessels ();
  int add_vessel (const int xc, const int yc, const VesselType tp);

  /* epithelial cells */
  int num_epi_in_initial_box;
  std::vector<Epithelial*> epithelium;
  int init_epithelium ();
  int add_epithelial (const int x, const int y, const EpithelialType tp, const EpithelialState st);
  void remove_epithelial (const int del_i);

  /* virions */
  int init_virions ();

  /* agents */
  std::vector<Agent*> all_agents;
  int init_all_agents ();

  /* convert number densities to particle numbers */
  double density_to_particle (Epithelial* epi, const bool infectious);

  /* templated agent routines */
  template <class T, typename T_TYP, typename T_ST>
  int init_agent (const bool fixed, const std::vector<std::pair<int,int>> fixed_coords, const int initial_number, const T_TYP tp, const T_ST st);

  template <class T, typename T_TYP, typename T_ST>
  int add_agent (const int x, const int y, const T_TYP tp, const T_ST st, const bool is_converted_agent, const bool use_output);

  /* non-templated agent routines */
  int remove_agent (const int del_i);
  int move_agent (const int agent_index, const int epithelial_index);

  /* compute viral spread */
  int compute_spread (const double time_n, const double time);

  /* compute epithelial death */
  int compute_epithelial_death (const double time);

  /* compute the action of the agents */
  int compute_agent_action (const double time);

  /* compute agent migration */
  int compute_agent_migration (const double time);

  /* compute agent recruitment */
  int compute_agent_recruitment ();

  /* templated routine for the recruitment of agents */
  template <class T, typename T_TYP, typename T_ST>
  int recruit_agents (Vessel* ves, const T_TYP tp, const T_ST st);

  /* find an available space */
  int find_space (Agent* agent_to_move, std::vector<std::pair<int,int>> neighbours);

  /* check a space */
  bool check_space (Agent* agent_to_move, std::pair<int,int> xy);

  /* templated output routine for agent codes */
  template <class T, typename T_TYP, typename T_ST>
  int get_agent_codes (int** grid_codes, const T_TYP tp, const T_ST st, const double time);

public:
  Environment (int* status);
  ~Environment ();

  /* advance the time step */
  int advance (const double time_n, const double time);

  /* write the consts */
  void write_consts (struct cam_output_files* files);

  /* write the data */
  void write_data (struct cam_output_files* files, const double time);
};

#endif /* __CAM_ENVIRONMENT_H__ */
