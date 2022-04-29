#include <cstdio>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <yaml-cpp/yaml.h>
#include <cmath>

#include "CAModel/cam_progress.h"
#include "CAModel/cam_error.h"
#include "CAModel/cam_output.h"
#include "CAModel/cam_parameters.h"
#include "CAModel/cam_environment.h"


/* output folder */
std::string output_folder;


/* CA model usage string, printed if bad args */
const char cam_usage_string[] =
	"camodel.exe <YAML parameter file> <output directory>\n";


/* function to check the correct YAML file extension */
int check_yaml_fextn (const char* fname);


int
main (int argc, char** argv)
{
	int status = 0;

  printf("\n===================================");
	printf("\nSARS-CoV-2 Cellular Automaton Model");
	printf("\n===================================");
	printf("\n\n");
	fflush(stdout);

  if (argc != 3) {
		printf("USAGE: \n\t %s\n\n", cam_usage_string);
		CAM_ERROR("incorrect number of arguments", 9);
  }

  /* check extension of given file name */
	status = check_yaml_fextn(argv[1]);
	CAM_ERROR_CHECK(status, CAM_SUCCESS);

	/* load the YAML file */
  YAML::Node cfg = YAML::LoadFile(argv[1]);

  /* set the output folder from input args */
  output_folder = argv[2];
	if (!strcmp(&(output_folder.back()), "/")) {
		/* remove trailing slash if it exists */
		output_folder.pop_back();
	}

	/* check that the output folder exists */
  struct stat info;
  if (stat(output_folder.c_str(), &info) != 0) {
		CAM_ERROR("output folder does not exist", 9);
  }

  /* open the output files */
  struct cam_output_files files;
  status = cam_output_open(output_folder, &files);
	CAM_ERROR_CHECK(status, CAM_SUCCESS);


  /*
   * =========================================================================
   * SETUP
   * =========================================================================
   */

	/* set the global parameters */
	status = cam_parameters_set_global(&cfg);
	if (status != CAM_SUCCESS) {
		fprintf(stderr, "\n");

		/* clear up */
		cam_output_close(&files);

		return status;
	}

	/* set the tissue environment parameters */
	status = cam_env_set_params(&cfg);
	if (status != CAM_SUCCESS) {
		fprintf(stderr, "\n");

		/* clear up */
		cam_output_close(&files);

		return status;
	}

	/* initialise tissue environment */
  Environment* t_env = new Environment(&status);
	if (status != CAM_SUCCESS) {
		fprintf(stderr, "\n");

		/* clear up */
		cam_output_close(&files);
		if (t_env != nullptr) {
			delete t_env;
		}

		return status;
	}


  /*
   * =========================================================================
   * SIMULATION
   * =========================================================================
   */

  printf("Time limit: %g [hr(s)]\n", time_limit);
  printf("Timestep: %g [hr(s)]\n\n", timestep);

	double time_n = 0.0;
  double time = 0.0;

  /* write the initial output */
  t_env->write_data(&files, time);

	/* start the clock */
  clock_t simulation_start_time = clock();

	unsigned int n = 1;
  while (time < time_limit) {
    /* Increment the time */
		time_n = time;
    time = ((double) n)*timestep;

		status = t_env->advance(time_n, time);
		if (status != CAM_SUCCESS) {
			fprintf(stderr, "\n");

			/* clear up */
			cam_output_close(&files);
			if (t_env != nullptr) {
				delete t_env;
			}

			return 9;
		}

    /* Write to files on the output interval */
    if ((int) round(time/timestep) % (int) round(output_interval/timestep) == 0) {
      PROGRESS(time);
      t_env->write_data(&files, time);
    }

		n++;
  }
	printf("\n");


  /* Output the global constants and close the output files */
	t_env->write_consts(&files);

	/* clear up */
	cam_output_close(&files);
	if (t_env != nullptr) {
		delete t_env;
	}

	/* Simulation complete - output time taken */
  double time_taken = (clock() - simulation_start_time)/CLOCKS_PER_SEC;
  printf("\nSIMULATION COMPLETE \n");
  printf("Time taken: %.2f [sec(s)]\tor %.2f [min(s)]\tor %.2f [hr(s)] \n", time_taken, time_taken/60.0, time_taken/60.0/60.0);

  return 0;
}


/*
 * =============================================================================
 * Check the YAML file extension
 * =============================================================================
 */
int
check_yaml_fextn (const char* fname)
{
	/* find last occurence of '.' in fname */
	const char* dot = strrchr(fname, '.');
	if (!dot || !strcmp(dot, fname)) {
		CAM_ERROR("no file extension", 9);
	}

	/* move the pointer forward to remove the '.' */
	dot++;

	/* known file extensions */
	const int yaml_n = 2;
	const char* yaml_extn[2] = {"yml", "yaml"};

	/* check file extension */
	for (int i = 0; i < yaml_n; i++) {
		const char* extn = yaml_extn[i];
		if (!strcmp(dot, extn))
			return CAM_SUCCESS;
	}

	/* unknown file extension */
	CAM_ERROR("unknown file extension", 9);
}
