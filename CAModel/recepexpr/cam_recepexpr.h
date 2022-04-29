#ifndef __CAM_RECEPEXPR_H__
#define __CAM_RECEPEXPR_H__

#include <yaml-cpp/yaml.h>

#include "../cam_error.h"


/* receptor expression type */
typedef struct {
  /* name of type */
  const char* name;

  /* dimension of system */
  unsigned int dim;

  /* functions */
  int (*set_params) (YAML::Node* cfg);                        // set the type specific parameters
  void* (*alloc_params) ();                                          // allocate the type specific data
  void* (*alloc_data) ();                                          // allocate the type specific data
  int (*init) (void* vpar, void* vdat, size_t dim);                       // initialise the type
  void (*input) (void* vpar, void* vinput, void* vifn);                   // set the input for the ODE solve
  int (*solve) (void* vpar, void* vdat, double t_start, double t_end, void* vupt, void* vsrc);    // solve the type
  void (*free) (void* vdat);                                  // free the type and type specific data
  void (*output) (void* vpar, void* vdat, const double time, FILE* outf);
  double (*get_sol) (void* vdat, const int index);

} cam_re_type;

/* available receptor expression types */
extern const cam_re_type* cam_re_type_pcv2_discon;


/* Set the receptor expression specific parameters */
int cam_recepexpr_set_params (YAML::Node* cfg, const cam_re_type* tp);


/*
 * =============================================================================
 * Receptor Expression
 * =============================================================================
 */
class RecepExpr
{
private:
  const cam_re_type* type;
  void* params; // type specific parameters
  void* data; // type specific data (hence void*)

  FILE* output;

public:
  RecepExpr (const cam_re_type* tp, bool use_output, int* status);
  ~RecepExpr ();

  void set_input (void* vin, void* vifn);
  int solve_system (double t_start, double t_end, void* uptake, void* source);

  void write_output (const double time);

  double get_data ();
};


/*
 * =============================================================================
 * Inline functions
 * =============================================================================
 */

/* set the input */
inline void
RecepExpr::set_input (void* vin, void* vifn)
{
  type->input(params, vin, vifn);
}


/* solve the ODE system */
inline int
RecepExpr::solve_system (double t_start, double t_end, void* uptake, void* source)
{
  int status = type->solve(params, data, t_start, t_end, uptake, source);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  return CAM_SUCCESS;
}


/* write the output */
inline void
RecepExpr::write_output (const double time)
{
  if (output == nullptr) {
    return;
  }

  type->output(params, data, time, output);
}

/* get the data */
inline double
RecepExpr::get_data ()
{
  /* internal bound receptors are at index 2 */
  double re_ib = type->get_sol(data, 2);

  return re_ib;
}

#endif /* __CAM_RECEPEXPR_H__ */
