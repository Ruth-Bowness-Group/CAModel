#ifndef __CAM_AGENT_H__
#define __CAM_AGENT_H__

#include <vector>
#include <functional>

#include "../reactdiff/cam_reactdiff.h"
#include "../cam_error.h"
#include "../recepexpr/cam_recepexpr.h"


/* enum of Agent Type */
enum class AgentType
{
  null = -1,
  virion = 100,
  macrophage = 200,
  neutrophil = 300,
  monocyte = 400,
  nk = 500
};


class Agent
{
private:
  /* coordinates */
  int x;
  int y;

  /* age */
  double age;
  double lifespan;
  double max_lifespan;

  /* area */
  double area;

  /* Agent type */
  AgentType ag_type;

protected:
  /* convert agent type */
  bool flag_conversion;

  /* bursting */
  bool flag_bursting;

  /* apoptosis */
  bool is_apoptosing_epithelial;
  Agent* apoptosing_agent;

  /* phagocytosis */
  bool is_phagocytosing_epithelial;

public:
  Agent (const int x, const int y, const AgentType tp, struct cam_rd_sys* rds);
  virtual ~Agent ();

  /* index */
  int index;

  /* apoptosis */
  double apoptosis_timescale;
  double apoptosis_start_time;
  bool is_being_apoptosed;
  int num_cells_apoptosing;
  virtual bool is_inflammatory ();
  virtual bool is_active ();
  virtual void start_apoptosis (const double time);
  virtual void end_apoptosis ();
  virtual void set_apoptosis_timescale (const double tms);

  /* phagocytosis */
  Agent* phagocytosing_agent;
  double phagocytosis_timescale;
  double phagocytosis_start_time;
  bool is_being_phagocytosed;
  Agent* phagocytosed_by;
  virtual bool is_apoptotic ();
  virtual bool is_phagocyte ();
  virtual void start_phagocytosis (const double time);
  virtual void end_phagocytosis ();
  virtual void set_phagocytosis_timescale (const double tms);

  /* reaction diffusion system */
  struct cam_rd_sys* rd_sys;

  /* coordinates */
  std::pair<int,int> get_coordinates ();
  void set_coordinates (const int cx, const int cy);

  /* age of cell in state */
  double get_age ();
  void set_age (const double amount);
  void increase_age (const double amount);

  /* get area */
  void set_area (const double area_to_set);
  double get_area ();

  /* agent type */
  AgentType get_agent_type ();
  void set_agent_type (const AgentType tp);
  bool check_agent_type (const AgentType tp);

  /* change state */
  virtual void change_state ();

  /* lifespan */
  double get_lifespan (const bool max);
  void set_lifespan (const double span, const bool max);
  bool exceeded_lifespan ();

  /* output code */
  virtual int get_code ();

  /* resident on cell */
  virtual void set_as_resident_on_cell (void* cell);
  virtual void* get_resident_on_cell ();

  virtual void dec_rds_values ();
  virtual void inc_rds_values ();

  virtual void action (const double time);
  virtual int migrate (const double dt, const double tm, std::function<int(Agent*, std::vector<std::pair<int,int>>)> fndsp, std::function<bool(Agent*, std::pair<int,int>)> chksp);
  virtual bool recruit (std::vector<std::pair<int,int>> nbrs);

  bool check_flag_conversion ();
  void switch_flag_conversion ();

  bool check_flag_bursting ();

  virtual double get_internalised ();
};


/*
 * =============================================================================
 * Inline functions
 * =============================================================================
 */

/* get the coordinates */
inline std::pair<int,int>
Agent::get_coordinates ()
{
  return std::make_pair(x, y);
}

/* set the coordinates */
inline void
Agent::set_coordinates (const int cx, const int cy)
{
  x = cx;
  y = cy;
}

/* get the age of the agent in current state */
inline double
Agent::get_age ()
{
  return age;
}

/* set the age of the agent */
inline void
Agent::set_age (const double amount)
{
  age = amount;
}

/* increase the age of the agent in current state */
inline void
Agent::increase_age (const double amount)
{
  age += amount;
}

/* set area */
inline void
Agent::set_area (const double area_to_set)
{
  area = area_to_set;
}

/* get area */
inline double
Agent::get_area ()
{
  return area;
}

/* get the agent type */
inline AgentType
Agent::get_agent_type ()
{
  return ag_type;
}

/* set the agent type */
inline void
Agent::set_agent_type (const AgentType tp)
{
  ag_type = tp;
}

/* check the agent type */
inline bool
Agent::check_agent_type (const AgentType tp)
{
  return ( (ag_type == tp) ? true : false );
}

/* change the state */
inline void
Agent::change_state ()
{
  /* nothing to do */
}

/* get the lifespan */
inline double
Agent::get_lifespan (const bool max)
{
  return ( (max) ? max_lifespan : lifespan );
}

/* set the lifespan */
inline void
Agent::set_lifespan (const double span, const bool max)
{
  lifespan = ( (span < CAM_TOL) ? 0.0 : span );
  if (max) {
    max_lifespan = lifespan;
  }
}

/* check whether exceeded lifespan */
inline bool
Agent::exceeded_lifespan ()
{
  return ( (age >= lifespan) ? true : false );
}

/* get the output code */
inline int
Agent::get_code ()
{
  return static_cast<int>(ag_type);
}

/* resident on cell */
inline void
Agent::set_as_resident_on_cell (void* cell)
{
  (void) cell; /* avoid unused parameter warning */
}

/* get the resident this agent is resident on */
inline void*
Agent::get_resident_on_cell ()
{
  return nullptr;
}

inline void
Agent::dec_rds_values ()
{
  // nothing to do
}

inline void
Agent::inc_rds_values ()
{
  // nothing to do
}

inline void
Agent::action (const double time)
{
  (void) time; /* avoid unused parameter warning */
}

// inline int
inline int
Agent::migrate (const double dt, const double tm, std::function<int(Agent*, std::vector<std::pair<int,int>>)> fndsp, std::function<bool(Agent*, std::pair<int,int>)> chksp)
{
  (void) dt; /* avoid unused parameter warning */
  (void) tm; /* avoid unused parameter warning */
  (void) fndsp; /* avoid unused parameter warning */
  (void) chksp; /* avoid unused parameter warning */

  return CAM_ERROR_MISC;
}

inline bool
Agent::recruit (std::vector<std::pair<int,int>> nbrs)
{
  (void) nbrs; /* avoid unused parameter warning */

  return false;
}

inline bool
Agent::check_flag_conversion ()
{
  return flag_conversion;
}

/* switch the flag for conversion */
inline void
Agent::switch_flag_conversion ()
{
  flag_conversion = !flag_conversion;
}

/* check flag bursting */
inline bool
Agent::check_flag_bursting ()
{
  return flag_bursting;
}

/* get the internalised species */
inline double
Agent::get_internalised ()
{
  /* nothing to do */

  CAM_WARN("agent has no get_internalised routines", 0.0);
}


/*
 * -----------------------------------------------------------------------------
 * Apoptosis routines
 * -----------------------------------------------------------------------------
 */

/* check whether agent is inflammatory */
inline bool
Agent::is_inflammatory ()
{
  return false;
}

/* check whether agent is active */
inline bool
Agent::is_active ()
{
  return false;
}

inline void
Agent::start_apoptosis (const double time)
{
  (void) time; /* avoid unused parameter warning */

  CAM_WARN_VOID("this agent has no start_apoptosis routine");
}

inline void
Agent::end_apoptosis ()
{
  /* nothing to do */
  CAM_WARN_VOID("this agent has no end_apoptosis routine");
}

inline void
Agent::set_apoptosis_timescale (const double tms)
{
  (void) tms; /* avoid unused parameter warning */

  CAM_WARN_VOID("this agent has no set_apoptosis_timescale routine");
}


/*
 * -----------------------------------------------------------------------------
 * Phagocytosis routines
 * -----------------------------------------------------------------------------
 */

/* check whether agent is apoptotic */
inline bool
Agent::is_apoptotic ()
{
  return false;
}

inline bool
Agent::is_phagocyte ()
{
  return false;
}

inline void
Agent::start_phagocytosis (const double time)
{
  (void) time; /* avoid unused parameter warning */
}

inline void
Agent::end_phagocytosis ()
{
  /* nothing to do */
}

inline void
Agent::set_phagocytosis_timescale (const double tms)
{
  (void) tms; /* avoid unused parameter warning */
}

#endif /* __CAM_AGENT_H__ */
