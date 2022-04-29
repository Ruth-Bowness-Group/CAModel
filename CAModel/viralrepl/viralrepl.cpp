#include "cam_viralrepl.h"
#include "../cam_error.h"
#include "../cam_parameters.h"


/*
 * =============================================================================
 * Set the viral replication specific parameters
 * =============================================================================
 */
int
cam_viralrepl_set_params (YAML::Node* cfg, const cam_vr_type* tp)
{
  if (tp == nullptr) {
    return CAM_SUCCESS;
  }

  /* call the set parameters specific to the type */
  int status = tp->set_params(cfg);
  CAM_ERROR_CHECK(status, CAM_SUCCESS);

  return CAM_SUCCESS;
}


/*
 * =============================================================================
 * ViralRepl constructor
 * =============================================================================
 */
ViralRepl::ViralRepl (const cam_vr_type* tp, bool use_output, int* status)
{
  if (tp == nullptr) {
    type = nullptr;
    data = nullptr;

    (*status) = CAM_SUCCESS;
    return;
  }

  /* set the viral replication type */
  type = tp;

  /* allocate the type specific params */
  params = type->alloc_params();
  if (params == nullptr) {
    (*status) = CAM_ERROR_VALUE;
    CAM_ERROR_VOID("cannot allocate type specific params");
  }

  /* allocate the type specific data */
  data = type->alloc_data();
  if (data == nullptr) {
    (*status) = CAM_ERROR_VALUE;
    CAM_ERROR_VOID("cannot allocate type specific data");
  }

  /* initialise the type */
  (*status) = type->init(params, data, type->dim);
  CAM_ERROR_CHECK_VOID((*status), CAM_SUCCESS);

  if (use_output) {
    const std::string outf = output_folder + std::string("/vr.txt");
    output = fopen(outf.c_str(), "w");
  } else {
    output = nullptr;
  }

  /* successful */
  (*status) = CAM_SUCCESS;
}


/*
 * =============================================================================
 * ViralRepl destructor
 * =============================================================================
 */
ViralRepl::~ViralRepl ()
{
  /* free the type specific data */
  if (type != nullptr) {
    type->free(data);
  }

  if (params != nullptr) {
    free(params);
  }

  if (data != nullptr) {
    free(data);
  }

  if (output != nullptr) {
    fclose(output);
  }
}
