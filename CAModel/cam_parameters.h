#ifndef __CAM_PARAMETERS_H__
#define __CAM_PARAMETERS_H__

#include <yaml-cpp/yaml.h>
#include <string>

#ifndef POS_ZERO
# define POS_ZERO (+0.0)
#endif

#ifndef NEG_ZERO
# define NEG_ZERO (-0.0)
#endif

#ifndef CAM_EPS
# define CAM_EPS (1e-15)
#endif

#ifndef CAM_TOL
# define CAM_TOL (1e-12)
#endif

#ifndef CAM_PI
# define CAM_PI (3.14159265358979323846264338328)
#endif

#ifndef CAM_CMP
# define CAM_CMP(x,y,e) ( ((y)-(x)) > (e) ) /* x < y */
#endif


/* output parameters */
extern std::string output_folder;
extern double output_interval;

/* simulation parameters */
extern double time_limit;
extern double timestep;

/* environment parameters */
const int grid_size = 200;

/* openmp threads - NOT USED */
extern int cam_omp_num_threads;

/* boolean flag for chemotactic migration and recruitment validation only */
const bool validation_chemotaxis_and_recruitment = false;


/* set the global parameters */
int cam_parameters_set_global (YAML::Node* cfg);

#endif /* __CAM_PARAMETERS_H__ */
