#ifndef __CAM_VIRALREPL_H__
#define __CAM_VIRALREPL_H__

#include <yaml-cpp/yaml.h>

#include "../cam_error.h"


/* viral replication type */
typedef struct {
  /* name of type */
  const char* name;

  /* dimension of system */
  unsigned int dim;

  /* functions */
  int (*set_params) (YAML::Node* cfg);                        // set the type specific parameters
  void* (*alloc_params) ();                                          // allocate the type specific params
  void* (*alloc_data) ();                                          // allocate the type specific data
  int (*init) (void* vpar, void* vdat, size_t dim);                       // initialise the type
  void (*input) (void* vpar, void* vinput, void* vifn);                   // set the input for the ODE solve
  int (*solve) (void* vpar, void* vdat, double t_start, double t_end, void* vexp);    // solve the type
  void (*free) (void* vdat);                                  // free the type and type specific data
  void (*output) (void* vpar, void* vdat, const double time, FILE* outf);
  double (*get_sol) (void* vdat, const int index);
  bool (*get_exporting) (void* vpar);

} cam_vr_type;

/* available viral repl types */
extern const cam_vr_type* cam_vr_type_pcv4_discon;


/* Set the viral replication specific parameters */
int cam_viralrepl_set_params (YAML::Node* cfg, const cam_vr_type* tp);


/*
 * =============================================================================
 * Viral Replication
 * =============================================================================
 */
class ViralRepl
{
private:
  const cam_vr_type* type;
  void* params; // type specific params (hence void*)
  void* data; // type specific data (hence void*)

  FILE* output;

public:
  ViralRepl (const cam_vr_type* tp, bool use_output, int* status);
  ~ViralRepl ();

  void set_input (void* vin, void* vifn);
  int solve_system (double t_start, double t_end, void* source);

  void write_output (const double time);

  bool get_param ();
  double get_data ();
};


/*
 * =============================================================================
 * Inline functions
 * =============================================================================
 */

/* set the input */
inline void
ViralRepl::set_input (void* vin, void* vifn)
{
  type->input(params, vin, vifn);
}


/* solve the ODE system */
inline int
ViralRepl::solve_system (double t_start, double t_end, void* source)
{
  int status = type->solve(params, data, t_start, t_end, source);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  return CAM_SUCCESS;
}

/* write the output */
inline void
ViralRepl::write_output (const double time)
{
  if (output == nullptr) {
    return;
  }

  type->output(params, data, time, output);
}

/* get the parameter */
inline bool
ViralRepl::get_param ()
{
  /* are we exporting virions */
  bool exporting = type->get_exporting(params);

  return exporting;
}

/* get the data */
inline double
ViralRepl::get_data ()
{
  /* assembled virions are at index 4 */
  double vr_a = type->get_sol(data, 4);

  return vr_a;
}

#endif /* __CAM_VIRALREPL_H__ */
