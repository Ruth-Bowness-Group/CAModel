#ifndef __CAM_EPITHELIAL_H__
#define __CAM_EPITHELIAL_H__

#include <yaml-cpp/yaml.h>

#include "../reactdiff/cam_reactdiff.h"
#include "../viralrepl/cam_viralrepl.h"
#include "../recepexpr/cam_recepexpr.h"
#include "../agent/cam_agent.h"


/* Set the Epithelial specific parameters */
int cam_epithelial_set_params (YAML::Node* cfg);


/* enum of Epithelial Type */
enum class EpithelialType
{
  epithelial = 10
};


/* enum of Epithelial State */
enum class EpithelialState
{
  healthy = 0,
  eclipse = 1,
  infected = 2,
  apoptotic = 3,
  phagocytosed = 4,
  burst = 5
};


class Epithelial
{
private:
  /* coordinates */
  int x;
  int y;

  /* age */
  double age;

  /* area */
  double area;
  double max_area;

  /* epithelial type */
  EpithelialType type;

  /* epithelial state */
  EpithelialState state;

  /* reaction diffusion system */
  struct cam_rd_sys* rd_sys;

  /* T1 IFN */
  double basal_secreted_ifn_1;
  double secreted_ifn_1;
  double uptake_ifn_1;
  double epithelial_infected_delay_ifn_1;

  /* IL-6 */
  double basal_secreted_il_6;
  double secreted_il_6;
  double uptake_il_6;

  /* IL-10 */
  double basal_secreted_il_10;
  double secreted_il_10;
  double uptake_il_10;

  /* IL-12 */
  double basal_secreted_il_12;
  double secreted_il_12;
  double uptake_il_12;

  /* G-CSF */
  double basal_secreted_g_csf;
  double secreted_g_csf;
  double uptake_g_csf;

  /* M-CSF */
  double basal_secreted_m_csf;
  double secreted_m_csf;
  double uptake_m_csf;

  /* Infected epithelial chemokine */
  double basal_secreted_epi_inf_ck;
  double secreted_epi_inf_ck;
  double uptake_epi_inf_ck;

  /* apoptosis chemokine */
  double secreted_apop_ck;
  double uptake_apop_ck;

  /* phagocytosis chemokine */
  double secreted_phago_ck;
  double uptake_phago_ck;

  /* intracellular viral replication */
  ViralRepl* vr;

  /* receptor expression */
  RecepExpr* re;

  /* resident cells */
  std::vector<Agent*> residents;

  /* flag to say whether cell is set to burst or apoptose */
  bool flag_bursting;
  bool flag_apoptosis;
  bool flag_phagocytosis;

  /* set the cytokine secretion */
  void set_secretion (const EpithelialState st);
  void regulate_secretion ();
  void update_secretions ();

  /* set the uptake */
  void set_uptake (const EpithelialState st);
  void regulate_uptake ();

public:
  Epithelial (const int x, const int y, const EpithelialType tp, const EpithelialState st, struct cam_rd_sys* rds, int* status);
  ~Epithelial ();

  /* index */
  int index; // set externally

  /* epithelial is being phagocytosed */
  bool is_being_phagocytosed;
  double phagocytosis_timescale;
  double phagocytosis_start_time;

  /* epithelial is being apoptosed */
  bool is_being_apoptosed;
  double apoptosis_timescale;
  double apoptosis_start_time;
  int num_cells_apoptosing;

  /* neighbours */
  std::vector<std::pair<int,int>> immediate_moore;
  std::vector<std::pair<int,int>> neighbours;

  /* coordinates */
  std::pair<int,int> get_coordinates ();
  void set_coordinates (const int cx, const int cy);

  /* age of cell in state */
  double get_age ();
  void set_age (const double amount);
  void increase_age (const double amount);

  /* get area */
  double get_max_area ();
  bool check_area (const double area_to_add);

  /* epithelial cell type */
  EpithelialType get_type ();
  void set_type (const EpithelialType tp);

  /* epithelial cell state */
  EpithelialState get_state ();
  void set_state (const EpithelialState st);
  bool check_state (const EpithelialState st);
  bool check_flag_bursting ();
  void switch_flag_bursting ();
  bool check_flag_apoptosis ();
  void switch_flag_apoptosis ();
  bool check_flag_phagocytosis ();
  void switch_flag_phagocytosis ();
  void change_state (const EpithelialState st);

  /* change the RDS values (sources/uptakes etc..) */
  void dec_rds_values ();
  void inc_rds_values ();

  /* output code */
  int get_code ();

  /* residents */
  int add_resident (Agent* res);
  std::vector<Agent*>* get_residents ();
  void get_residents (AgentType tp, std::vector<Agent*>* res);
  int remove_resident (Agent* res);
  bool has_residents (AgentType tp);
  bool has_residents_inflammatory ();
  bool has_residents_apoptotic ();

  /* apoptosis */
  void start_apoptosis (const double time);
  void end_apoptosis ();
  void set_apoptosis_timescale (const double tms);

  /* phagocytosis */
  void start_phagocytosis (const double time);
  void end_phagocytosis ();
  void set_phagocytosis_timescale (const double tms);

  /* replication */
  int replicate (double nvirs, double time_n, double time);

  /* get interior components */
  double get_interior ();

  /* write VR and RE output */
  void write_vr_output (const double time);
  void write_re_output (const double time);
};


/*
 * =============================================================================
 * Inline functions
 * =============================================================================
 */

/* get the coordinates */
inline std::pair<int,int>
Epithelial::get_coordinates ()
{
  return std::make_pair(x, y);
}

/* set the coordinates */
inline void
Epithelial::set_coordinates (const int cx, const int cy)
{
  x = cx;
  y = cy;
}

/* get the age of the epithelial cell in current state */
inline double
Epithelial::get_age ()
{
  return age;
}

/* set the age of the epithelial cell */
inline void
Epithelial::set_age (const double amount)
{
  age = amount;
}

/* increase the age of the epithelial cell in current state */
inline void
Epithelial::increase_age (const double amount)
{
  age += amount;
}

/* check whether agent can be added based on area */
inline bool
Epithelial::check_area (const double area_to_add)
{
  double tmp_area = area + area_to_add;
  return ( (CAM_CMP(tmp_area, max_area, CAM_TOL)) ? true : false ); /* area + area_to_add < max_area */
}

/* get the area */
inline double
Epithelial::get_max_area ()
{
  return max_area;
}

/* get epithelial type */
inline EpithelialType
Epithelial::get_type ()
{
  return type;
}

/* set the epithelial type */
inline void
Epithelial::set_type (const EpithelialType tp)
{
  type = tp;
}

/* get epithelial state */
inline EpithelialState
Epithelial::get_state ()
{
  return state;
}

/* check cell has burst */
inline bool
Epithelial::check_state (const EpithelialState st)
{
  return ( (state == st) ? true : false );
}

inline void
Epithelial::dec_rds_values ()
{
  /* coords */
  std::pair<int, int> xy = get_coordinates();

  /* need to decrease the cytokine secretions and uptakes */

  /* T1 IFN */
  rd_sys->ifn_1->dec_source(xy.first, xy.second, secreted_ifn_1, false);
  rd_sys->ifn_1->dec_uptake(xy.first, xy.second, uptake_ifn_1, true);

  /* IL-6 */
  rd_sys->il_6->dec_source(xy.first, xy.second, secreted_il_6, false);
  rd_sys->il_6->dec_uptake(xy.first, xy.second, uptake_il_6, true);

  /* IL-10 */
  rd_sys->il_10->dec_source(xy.first, xy.second, secreted_il_10, false);
  rd_sys->il_10->dec_uptake(xy.first, xy.second, uptake_il_10, true);

  /* IL-12 */
  rd_sys->il_12->dec_source(xy.first, xy.second, secreted_il_12, false);
  rd_sys->il_12->dec_uptake(xy.first, xy.second, uptake_il_12, true);

  /* G-CSF */
  rd_sys->g_csf->dec_source(xy.first, xy.second, secreted_g_csf, false);
  rd_sys->g_csf->dec_uptake(xy.first, xy.second, uptake_g_csf, true);

  /* M-CSF */
  rd_sys->m_csf->dec_source(xy.first, xy.second, secreted_m_csf, false);
  rd_sys->m_csf->dec_uptake(xy.first, xy.second, uptake_m_csf, true);

  /* epi_inf_ck */
  rd_sys->epi_inf_ck->dec_source(xy.first, xy.second, secreted_epi_inf_ck, false);
  rd_sys->epi_inf_ck->dec_uptake(xy.first, xy.second, uptake_epi_inf_ck, true);

  /* phago_ck */
  rd_sys->phago_ck->dec_source(xy.first, xy.second, secreted_phago_ck, false);
  rd_sys->phago_ck->dec_uptake(xy.first, xy.second, uptake_phago_ck, true);
}

inline void
Epithelial::inc_rds_values ()
{
  /* coords */
  std::pair<int, int> xy = get_coordinates();

  /* increase the secretions and uptakes */

  /* T1 IFN */
  rd_sys->ifn_1->inc_source(xy.first, xy.second, secreted_ifn_1, false);
  rd_sys->ifn_1->inc_uptake(xy.first, xy.second, uptake_ifn_1, true);

  /* IL-6 */
  rd_sys->il_6->inc_source(xy.first, xy.second, secreted_il_6, false);
  rd_sys->il_6->inc_uptake(xy.first, xy.second, uptake_il_6, true);

  /* IL-10 */
  rd_sys->il_10->inc_source(xy.first, xy.second, secreted_il_10, false);
  rd_sys->il_10->inc_uptake(xy.first, xy.second, uptake_il_10, true);

  /* IL-12 */
  rd_sys->il_12->inc_source(xy.first, xy.second, secreted_il_12, false);
  rd_sys->il_12->inc_uptake(xy.first, xy.second, uptake_il_12, true);

  /* G-CSF */
  rd_sys->g_csf->inc_source(xy.first, xy.second, secreted_g_csf, false);
  rd_sys->g_csf->inc_uptake(xy.first, xy.second, uptake_g_csf, true);

  /* M-CSF */
  rd_sys->m_csf->inc_source(xy.first, xy.second, secreted_m_csf, false);
  rd_sys->m_csf->inc_uptake(xy.first, xy.second, uptake_m_csf, true);

  /* epi_inf_ck */
  rd_sys->epi_inf_ck->inc_source(xy.first, xy.second, secreted_epi_inf_ck, false);
  rd_sys->epi_inf_ck->inc_uptake(xy.first, xy.second, uptake_epi_inf_ck, true);

  /* phago_ck */
  rd_sys->phago_ck->inc_source(xy.first, xy.second, secreted_phago_ck, false);
  rd_sys->phago_ck->inc_uptake(xy.first, xy.second, uptake_phago_ck, true);
}

/* check whether cell is flagged to burst */
inline bool
Epithelial::check_flag_bursting ()
{
  return flag_bursting;
}

/* switch the flag for bursting */
inline void
Epithelial::switch_flag_bursting ()
{
  flag_bursting = !flag_bursting;
}

/* check whether cell is flagged to apoptose */
inline bool
Epithelial::check_flag_apoptosis ()
{
  return flag_apoptosis;
}

/* switch the flag for apoptosis */
inline void
Epithelial::switch_flag_apoptosis ()
{
  flag_apoptosis = !flag_apoptosis;
}

/* check whether cell is flagged to phagocytose */
inline bool
Epithelial::check_flag_phagocytosis ()
{
  return flag_phagocytosis;
}

/* switch the flag for phagocytosis */
inline void
Epithelial::switch_flag_phagocytosis ()
{
  flag_phagocytosis = !flag_phagocytosis;
}

/* get the output grid code */
inline int
Epithelial::get_code ()
{
  return (static_cast<int>(type) + static_cast<int>(state));
}

/* get the residents */
inline std::vector<Agent*>*
Epithelial::get_residents ()
{
  return &residents;
}

inline void
Epithelial::write_re_output (const double time)
{
  re->write_output(time);
}

inline void
Epithelial::write_vr_output (const double time)
{
  vr->write_output(time);
}

#endif /* __CAM_EPITHELIAL_H__ */
