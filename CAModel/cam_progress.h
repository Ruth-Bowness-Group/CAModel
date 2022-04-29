#ifndef __CAM_PROGRESS_H__
#define __CAM_PROGRESS_H__

#include <stdio.h>

#define PROGRESS(tm) \
  do { \
    fprintf(stdout, "Time: %g [hr(s)]\n", (tm)); \
    fflush(stdout); \
  } while(0)


#endif /* __CAM_PROGRESS_H__ */
