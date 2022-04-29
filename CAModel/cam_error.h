#ifndef __CAM_ERROR_H__
#define __CAM_ERROR_H__

/*******************************************************************************
 * Error codes
 ******************************************************************************/

/* enum for warning (negative) and error (positive) codes */
enum {
  CAM_WARN_MISC = -3,                /* misc warning */
  CAM_WARN_AREA = -2,                /* area */
  CAM_WARN_EXISTS = -1,              /* exists */
  CAM_SUCCESS = 0,                   /* success */
  CAM_ERROR_IO = 1,                  /* failure */
  CAM_ERROR_NOMEM = 2,               /* allocation */
  CAM_ERROR_BADARGS = 3,             /* bad args to a routine */
  CAM_ERROR_VALUE = 4,               /* value */
  CAM_ERROR_EXISTS = 5,              /* exists */
  CAM_ERROR_GSL = 6,                 /* gnu scientific library */
  CAM_ERROR_MISC = 9                 /* misc error not matching above */
};


/*******************************************************************************
 * Macros to handle the error
 ******************************************************************************/

/* handle the error and return an error code */
#define CAM_ERROR(em, ec) \
  do { \
    fprintf(stderr, "cam: %s:%s:%d => ERROR: %s\n", __FILE__, __func__, __LINE__, (em)); \
    return (ec); \
  } while (0)

/* handle the warning and return a warning code */
#define CAM_WARN(em, ec) \
  do { \
    fprintf(stderr, "cam: %s:%s:%d => WARN: %s\n", __FILE__, __func__, __LINE__, (em)); \
    return (ec); \
  } while (0)

#define CAM_WARN_VOID(em) \
  do { \
    fprintf(stderr, "cam: %s:%s:%d => WARN: %s\n", __FILE__, __func__, __LINE__, (em)); \
    return; \
  } while (0)

/* handle the error without return */
#define CAM_ERROR_NORETURN(em) \
  fprintf(stderr, "cam: %s:%s:%d => ERROR: %s\n", __FILE__, __func__, __LINE__, (em));

/* handle the error and return from a void function */
#define CAM_ERROR_VOID(em) \
  do { \
    fprintf(stderr, "cam: %s:%s:%d => ERROR: %s\n", __FILE__, __func__, __LINE__, (em)); \
    return; \
  } while (0)

/* propagate the error to calling routines */
#define CAM_ERROR_TRACE(ec) \
  do { \
    fprintf(stderr, "cam: %s:%s:%d\n", __FILE__, __func__, __LINE__); \
    return (ec); \
  } while (0)

/* propagate the error to calling routines */
#define CAM_ERROR_TRACE_VOID(ec) \
  do { \
    fprintf(stderr, "cam: %s:%s:%d\n", __FILE__, __func__, __LINE__); \
    return; \
  } while (0)

/* check return status of a routine; if fail trace error */
#define CAM_ERROR_CHECK(stat, ec) \
  do { \
    if ((stat) != (ec)) \
      CAM_ERROR_TRACE(stat); \
  } while (0)

/* check return status of a routine; if fail trace error */
#define CAM_ERROR_CHECK_VOID(stat, ec) \
  do { \
    if ((stat) != (ec)) \
      CAM_ERROR_TRACE_VOID(stat); \
  } while (0)

/* check output file opened; if fail trace error */
#define CAM_ERROR_CHECK_IO(fp, em) \
  do { \
    if ((fp) == nullptr) { \
      fprintf(stderr, "cam: %s:%s:%d => ERROR: cannot open output file (%s)\n", __FILE__, __func__, __LINE__, (em)); \
      return (CAM_ERROR_IO); \
    } \
  } while (0)

#endif /* __CAM_ERROR_H__ */
