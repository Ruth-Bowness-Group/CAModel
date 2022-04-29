#include "cam_parameters.h"
#include "cam_error.h"


/* output parameters */
double output_interval;

/* simulation parameters */
double time_limit;
double timestep;

/* openmp threads - NOT USED */
int cam_omp_num_threads;


/*
 * =============================================================================
 * Set the global parameters
 * =============================================================================
 */
int
cam_parameters_set_global (YAML::Node* cfg)
{
	try {
		output_interval = (*cfg)["output_interval"].as<double>();
	  time_limit = (*cfg)["time_limit"].as<double>();
	  timestep = (*cfg)["timestep"].as<double>();
		cam_omp_num_threads = (*cfg)["omp_num_threads"].as<int>(); // NOT USED

	} catch (std::exception &e) {
		CAM_ERROR(e.what(), CAM_ERROR_IO);
	}

	return CAM_SUCCESS;
}
