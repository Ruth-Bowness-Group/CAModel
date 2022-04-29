#ifndef __CAM_REACTDIFF_H__
#define __CAM_REACTDIFF_H__

#include <yaml-cpp/yaml.h>

#include "../cam_parameters.h"


/* Set the reaction diffusion specific parameters */
int cam_rd_set_params (YAML::Node* cfg);


/* available reaction diffusion types */
enum class RDType
{
  ifn_1 = 0,
  il_6 = 1,
  il_10 = 2,
  il_12 = 3,
  g_csf = 4,
  m_csf = 5,
  epi_inf_ck = 6,
  apop_ck = 7,
  phago_ck = 8,
  virus_inf = 9,
  virus_ninf = 10
};


class ReactDiff
{
private:
  RDType type;

  double max_conc;
  double conc[grid_size][grid_size];
  double basal_source[grid_size][grid_size];
  double source[grid_size][grid_size];
  double burst_source[grid_size][grid_size];
  double diffusion[grid_size][grid_size];
  double basal_uptake[grid_size][grid_size];
  double uptake[grid_size][grid_size];
  double decay[grid_size][grid_size];

  /* Diffusion */
  double calc_diffusion (double v_c, double v_n, double v_e, double v_s, double v_w,
                         double d_c, double d_n, double d_e, double d_s, double d_w);

  double calc_diffusion_bdry (double v_c, double v_n, double v_e, double v_s, double v_w, double d_c);

public:
  ReactDiff (RDType tp, int* status);
  ~ReactDiff ();

  /* Type of ReactDiff system */
  void set_type (RDType tp);
  RDType get_type ();

  /* Get, set, increase the concentration at coords */
  double get_scaled (const int x, const int y);
  double get_conc (const int x, const int y);
  void inc_conc (const int x, const int y, const double c);
  void dec_conc (const int x, const int y, const double c);
  double get_max ();

  /* Get, set, increase, decrease source at coords */
  double get_source (const int x, const int y, const bool cst);
  double get_burst_source (const int x, const int y);
  void set_source (const int x, const int y, const double s, const bool cst);
  void set_burst_source (const int x, const int y, const double s);
  void inc_source (const int x, const int y, const double s, const bool cst);
  void inc_burst_source (const int x, const int y, const double s);
  void dec_source (const int x, const int y, const double s, const bool cst);
  void dec_burst_source (const int x, const int y, const double s);
  void reset (const int x, const int y);

  /* Get, set, increase, decrease uptake at coords */
  double get_uptake (const int x, const int y, const bool cst);
  void set_uptake (const int x, const int y, const double up, const bool cst);
  void inc_uptake (const int x, const int y, const double up, const bool cst);
  void dec_uptake (const int x, const int y, const double up, const bool cst);

  /* Get, set, increase, decrease decay at coords */
  double get_decay (const int x, const int y);
  void set_decay (const int x, const int y, const double d);
  void inc_decay (const int x, const int y, const double d);
  void dec_decay (const int x, const int y, const double d);

  /* set the initial conditions */
  void set_initial_conditions (const int x, const int y, const double val);

  /* Advance the reaction diffusion system */
  int advance ();
};


/*
 * =============================================================================
 * inline functions
 * =============================================================================
 */

/* set reaction diffusion type */
inline void
ReactDiff::set_type (RDType tp)
{
  type = tp;
}

/* get the reaction diffusion type */
inline RDType
ReactDiff::get_type ()
{
  return type;
}

/* get the scaled concentration */
inline double
ReactDiff::get_scaled (const int x, const int y)
{
  return ((max_conc > CAM_TOL) ? (conc[x][y] / max_conc) : -1.0);
}

/* get the concentration */
inline double
ReactDiff::get_conc (const int x, const int y)
{
  return conc[x][y];
}

/* reset the sources, uptakes and decays */
inline void
ReactDiff::reset (const int x, const int y)
{
  basal_source[x][y] = 0.0;
  source[x][y] = 0.0;
  burst_source[x][y] = 0.0;

  basal_uptake[x][y] = 0.0;
  uptake[x][y] = 0.0;
}

/* increase the concentration */
inline void
ReactDiff::inc_conc (const int x, const int y, const double c)
{
  conc[x][y] += c;
}

/* decrease the concentration */
inline void
ReactDiff::dec_conc (const int x, const int y, const double c)
{
  conc[x][y] -= c;
}

/* get the maximum conc */
inline double
ReactDiff::get_max ()
{
  return max_conc;
}

/* get the source */
inline double
ReactDiff::get_source (const int x, const int y, const bool cst)
{
  return ( (cst) ? source[x][y] : basal_source[x][y] );
}

/* get the burst source */
inline double
ReactDiff::get_burst_source (const int x, const int y)
{
  return burst_source[x][y];
}

/* set the source */
inline void
ReactDiff::set_source (const int x, const int y, const double s, const bool cst)
{
  if (cst) {
    source[x][y] = s;
    return;
  }

  basal_source[x][y] = s;
}

/* set the burst source */
inline void
ReactDiff::set_burst_source (const int x, const int y, const double s)
{
  burst_source[x][y] = s;
}

/* increase source at coords */
inline void
ReactDiff::inc_source (const int x, const int y, const double s, const bool cst)
{
  if (cst) {
    source[x][y] += s;
    return;
  }

  basal_source[x][y] += s;
}

/* increase the burst source */
inline void
ReactDiff::inc_burst_source (const int x, const int y, const double s)
{
  burst_source[x][y] += s;
}

/* decrease source at coords */
inline void
ReactDiff::dec_source (const int x, const int y, const double s, const bool cst)
{
  if (cst) {
    source[x][y] -= s;

    if (source[x][y] < CAM_TOL) {
      source[x][y] = 0.0;
    }

    return;
  }

  basal_source[x][y] -= s;

  if (basal_source[x][y] < CAM_TOL) {
    basal_source[x][y] = 0.0;
  }
}

/* decrease the burst source */
inline void
ReactDiff::dec_burst_source (const int x, const int y, const double s)
{
  burst_source[x][y] -= s;

  if (burst_source[x][y] < CAM_TOL) {
    burst_source[x][y] = 0.0;
  }
}

/* get the uptake */
inline double
ReactDiff::get_uptake (const int x, const int y, const bool cst)
{
  return ( (cst) ? uptake[x][y] : basal_uptake[x][y] );
}

/* set the uptake */
inline void
ReactDiff::set_uptake (const int x, const int y, const double up, const bool cst)
{
  if (cst) {
    uptake[x][y] = up;
    return;
  }

  basal_uptake[x][y] = up;
}

/* increase uptake at coords */
inline void
ReactDiff::inc_uptake (const int x, const int y, const double up, const bool cst)
{
  if (cst) {
    uptake[x][y] += up;
    return;
  }

  basal_uptake[x][y] += up;
}

/* decrease uptake at coords */
inline void
ReactDiff::dec_uptake (const int x, const int y, const double up, const bool cst)
{
  if (cst) {
    uptake[x][y] -= up;

    if (uptake[x][y] < CAM_TOL) {
      uptake[x][y] = 0.0;
    }

    return;
  }

  basal_uptake[x][y] -= up;

  if (basal_uptake[x][y] < CAM_TOL) {
    basal_uptake[x][y] = 0.0;
  }
}

/* get the decay */
inline double
ReactDiff::get_decay (const int x, const int y)
{
  return decay[x][y];
}

/* set the decay */
inline void
ReactDiff::set_decay (const int x, const int y, const double d)
{
  decay[x][y] = d;
}

/* increase decay at coords */
inline void
ReactDiff::inc_decay (const int x, const int y, const double d)
{
  decay[x][y] += d;
}

/* decrease decay at coords */
inline void
ReactDiff::dec_decay (const int x, const int y, const double d)
{
  decay[x][y] -= d;

  if (decay[x][y] < CAM_TOL) {
    decay[x][y] = 0.0;
  }
}

/* set the initial conditions */
inline void
ReactDiff::set_initial_conditions (const int x, const int y, const double val)
{
  /* set concentration */
  conc[x][y] = val;
}


/*
 * =============================================================================
 * available reaction diffusion equations
 * =============================================================================
 */
struct cam_rd_sys
{
  const int n_rd = 11;
  ReactDiff* ifn_1 = nullptr;
  ReactDiff* il_6 = nullptr;
  ReactDiff* il_10 = nullptr;
  ReactDiff* il_12 = nullptr;
  ReactDiff* g_csf = nullptr;
  ReactDiff* m_csf = nullptr;
  ReactDiff* epi_inf_ck = nullptr;
  ReactDiff* apop_ck = nullptr;
  ReactDiff* phago_ck = nullptr;
  ReactDiff* virus_inf = nullptr;
  ReactDiff* virus_ninf = nullptr;
};

#endif /* __CAM_REACTDIFF_H__ */
