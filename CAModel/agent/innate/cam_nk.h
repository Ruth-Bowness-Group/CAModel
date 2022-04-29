#ifndef __CAM_NK_H__
#define __CAM_NK_H__

#include <cstdio>
#include <yaml-cpp/yaml.h>
#include <functional>

#include "../cam_agent.h"
#include "../../reactdiff/cam_reactdiff.h"
#include "../../epithelial/cam_epithelial.h"

/* Set the nk specific parameters */
int cam_nk_set_params (YAML::Node* cfg);

/* enum of NK Type */
enum class NKType
{
  cytotoxic = 0
};


/* enum of NK State */
enum class NKState
{
  null = -1,
  resting = 0,
  active = 1,
  apoptotic = 2,
  phagocytosed = 3
};


/*
 * =============================================================================
 * Natural Killer Cells
 * =============================================================================
 */
class NK : public Agent
{
private:
  FILE* output;
  FILE* uptake_output;
  FILE* lifespan_output;
  FILE* mrate_output;

  /* nk type */
  NKType nk_type;

  /* nk state */
  NKState nk_state;

  /* nk is resident on cell */
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

  /* movement */
  double movement_rate;
  double movement_bias;

  /* set the cytokine secretion */
  void set_secretion (const NKState st);
  void regulate_secretion ();
  void update_secretions ();

  /* lifespan */
  void update_lifespan (const NKState st);

  /* set the uptake */
  void set_uptake (const NKState st);
  void regulate_uptake ();

  /* movement */
  void set_movement (const NKState st);
  void update_movement (const NKState st);

  /* apoptosis */
  int apoptose_epithelial (const double time);
  void apoptose_immune (const double time);

public:
    NK (const int x, const int y, const NKType tp, const NKState st,
        struct cam_rd_sys* rds, const bool use_output, int* status);
    ~NK ();

    /* nk type */
    NKType get_type ();
    void set_type (const NKType tp);

    /* nk cell state */
    NKState get_state ();
    void set_state (const NKState st);
    void change_state (const NKState st);

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

/* get the nk type */
inline NKType
NK::get_type ()
{
  return nk_type;
}

/* set the nk type */
inline void
NK::set_type (const NKType tp)
{
  nk_type = tp;
}

/* get the nk state */
inline NKState
NK::get_state ()
{
  return nk_state;
}

/* get the output code */
inline int
NK::get_code ()
{
  return (static_cast<int>(Agent::get_code()) + static_cast<int>(nk_type) + static_cast<int>(nk_state));
}

/* set the nk as resident on a cell */
inline void
NK::set_as_resident_on_cell (void* cell)
{
  resident_on_cell = (Epithelial*) cell;
}

/* get the cell the nk is resident on  */
inline void*
NK::get_resident_on_cell ()
{
  return resident_on_cell;
}

inline void
NK::dec_rds_values ()
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

  /* apop_ck */
  rd_sys->apop_ck->dec_source(xy.first, xy.second, secreted_apop_ck, false);
  rd_sys->apop_ck->dec_uptake(xy.first, xy.second, uptake_apop_ck, true);
}

inline void
NK::inc_rds_values ()
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

  /* apop_ck */
  rd_sys->apop_ck->inc_source(xy.first, xy.second, secreted_apop_ck, false);
  rd_sys->apop_ck->inc_uptake(xy.first, xy.second, uptake_apop_ck, true);
}


/*
 * -----------------------------------------------------------------------------
 * Apoptosis routines
 * -----------------------------------------------------------------------------
 */

/* check whether agent is inflammatory */
inline bool
NK::is_inflammatory ()
{
  return false;
}

inline bool
NK::is_active ()
{
  return ( (nk_state == NKState::active) ? true : false );
}

/* check whether agent is apoptotic */
inline bool
NK::is_apoptotic ()
{
  return ( (nk_state == NKState::apoptotic) ? true : false );
}

inline bool
NK::is_phagocyte ()
{
  return false;
}

inline void
NK::write_output (const double time)
{
  if (output != nullptr) {
    /* coords */
    std::pair<int, int> xy = get_coordinates();

    fprintf(output, "%d %d\n", xy.first, xy.second);
  }

  if (uptake_output != nullptr) {
    fprintf(uptake_output, "%.5e %g\n", time, 0.0);
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
NK::get_internalised ()
{
  CAM_WARN("NK cells have no get_internalised routine", 0.0);
}

#endif /* __CAM_NK_H__ */
