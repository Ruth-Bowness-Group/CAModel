#ifndef __CAM_MONOCYTE_H__
#define __CAM_MONOCYTE_H__

#include <cstdio>
#include <yaml-cpp/yaml.h>
#include <functional>

#include "../cam_agent.h"
#include "../../reactdiff/cam_reactdiff.h"
#include "../../epithelial/cam_epithelial.h"

/* Set the monocyte specific parameters */
int cam_monocyte_set_params (YAML::Node* cfg);

/* enum of Monocyte Type */
enum class MonocyteType
{
  null = 0
};


/* enum of Monocyte State */
enum class MonocyteState
{
  null = -1,
  resting = 0,
  active = 1,
  apoptotic = 2,
  phagocytosed = 3
};


/*
 * =============================================================================
 * Monocyte
 * =============================================================================
 */
class Monocyte : public Agent
{
private:
  FILE* output;
  FILE* uptake_output;
  FILE* lifespan_output;
  FILE* mrate_output;

  /* monocyte type */
  MonocyteType mono_type;

  /* monocyte state */
  MonocyteState mono_state;

  /* monocyte is resident on cell */
  Epithelial* resident_on_cell;

  /* T1 IFN */
  double basal_secreted_ifn_1;
  double secreted_ifn_1;
  double uptake_ifn_1;

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

  /* infectious virus */
  double uptake_virus_inf;
  double decay_virus_inf;
  double total_virus_uptaken;

  /* movement */
  double movement_rate;
  double movement_bias;

  /* set the cytokine secretion */
  void set_secretion (const MonocyteState st);
  void regulate_secretion ();
  void update_secretions ();

  /* lifespan */
  void update_lifespan (const MonocyteState st);

  /* set the uptake */
  void set_uptake (const MonocyteState st);
  void regulate_uptake ();

  /* set the movement */
  void set_movement (const MonocyteState st);
  void update_movement (const MonocyteState st);

  int phagocytose_epithelial (const double time);
  void phagocytose_immune (const double time);

public:
    Monocyte (const int x, const int y, const MonocyteType tp, const MonocyteState st,
              struct cam_rd_sys* rds, const bool use_output, int* status);
    ~Monocyte ();

    /* monocyte type */
    MonocyteType get_type ();
    void set_type (const MonocyteType tp);

    /* monocyte cell state */
    MonocyteState get_state ();
    void set_state (const MonocyteState st);
    void change_state (const MonocyteState st);

    /* output code */
    int get_code ();

    /* resident on cell */
    void set_as_resident_on_cell (void* cell);
    void* get_resident_on_cell ();

    /* change the RDS values (sources/uptakes etc..) */
    void dec_rds_values ();
    void inc_rds_values ();

    /* action */
    void action (const double time);

    /* migration */
    int migrate (const double dt, const double tm, std::function<int(Agent*, std::vector<std::pair<int,int>>)> fndsp, std::function<bool(Agent*, std::pair<int,int>)> chksp);

    /* recruitment */
    bool recruit (std::vector<std::pair<int,int>> nbrs);

    /* apoptosis */
    bool is_inflammatory ();
    bool is_active ();
    void start_apoptosis (const double time);
    void end_apoptosis ();
    void set_apoptosis_timescale (const double tms);

    /* phagocytosis */
    bool is_apoptotic ();
    bool is_phagocyte ();
    void start_phagocytosis (const double time);
    void end_phagocytosis ();
    void set_phagocytosis_timescale (const double tms);

    /* output */
    void write_output(const double time);

    /* internalised */
    double get_internalised ();
};


/*
 * =============================================================================
 * Inline functions
 * =============================================================================
 */

/* get the monocyte type */
inline MonocyteType
Monocyte::get_type ()
{
  return mono_type;
}

/* set the monocyte type */
inline void
Monocyte::set_type (const MonocyteType tp)
{
  mono_type = tp;
}

/* get the monocyte state */
inline MonocyteState
Monocyte::get_state ()
{
  return mono_state;
}

/* get the output code */
inline int
Monocyte::get_code ()
{
  return (static_cast<int>(Agent::get_code()) + static_cast<int>(mono_type) + static_cast<int>(mono_state));
}

/* set the monocyte as resident on a cell */
inline void
Monocyte::set_as_resident_on_cell (void* cell)
{
  resident_on_cell = (Epithelial*) cell;
}

/* get the cell the monocyte is resident on  */
inline void*
Monocyte::get_resident_on_cell ()
{
  return resident_on_cell;
}

inline void
Monocyte::dec_rds_values ()
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

  /* Virus */
  double epi_area = ( (resident_on_cell == nullptr) ? 1.0 : resident_on_cell->get_max_area());
  rd_sys->virus_inf->dec_uptake(xy.first, xy.second, uptake_virus_inf / epi_area, false);
}

inline void
Monocyte::inc_rds_values ()
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

  /* virus */
  double epi_area = ( (resident_on_cell == nullptr) ? 1.0 : resident_on_cell->get_max_area());
  rd_sys->virus_inf->inc_uptake(xy.first, xy.second, uptake_virus_inf / epi_area, false);
}


/*
 * -----------------------------------------------------------------------------
 * Apoptosis routines
 * -----------------------------------------------------------------------------
 */

/* check whether agent is inflammatory */
inline bool
Monocyte::is_inflammatory ()
{
  return false;
}

inline bool
Monocyte::is_active ()
{
  return ( (mono_state == MonocyteState::active) ? true : false );
}

/* check whether monocyte is apoptotic */
inline bool
Monocyte::is_apoptotic ()
{
  return ( (mono_state == MonocyteState::apoptotic) ? true : false );
}

inline bool
Monocyte::is_phagocyte ()
{
  return true;
}

inline void
Monocyte::write_output (const double time)
{
  if (output != nullptr) {
    /* coords */
    std::pair<int, int> xy = get_coordinates();

    fprintf(output, "%d %d\n", xy.first, xy.second);
  }

  if (uptake_output != nullptr) {
    fprintf(uptake_output, "%.5e %g\n", time, total_virus_uptaken);
  }

  if (lifespan_output != nullptr) {
    double tmp = get_lifespan(false);
    fprintf(lifespan_output, "%.5e %g\n", time, tmp);
  }

  if (mrate_output != nullptr) {
    fprintf(mrate_output, "%.5e %g\n", time, movement_rate);
  }
}

/* get internalised */
inline double
Monocyte::get_internalised ()
{
  return total_virus_uptaken;
}

#endif /* __CAM_MONOCYTE_H__ */
