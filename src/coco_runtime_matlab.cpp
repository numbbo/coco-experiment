/**
 * @file coco_runtime_matlab.cpp
 * @brief Specific COCO runtime implementation for the Matlab language
 * that replaces coco_runtime_c.c with the Matlab-specific counterparts.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cstdarg>

#include <mex.h>

#include "coco.h"
#include "coco_utilities.cpp"

void coco_error(char const* message, ...) {
  va_list args;

  char buffer[500];
  int n;

  mexPrintf("COCO FATAL ERROR: ");
  va_start(args, message);
  n = vsprintf(buffer, message, args);
  va_end(args);
  mexPrintf(buffer);

  mexPrintf("\n");
  exit(EXIT_FAILURE);
}

void coco_warning(char const* message, ...) {
  va_list args;

  char buffer[500];
  int n;

  if (coco_log_level >= COCO_WARNING) {
    mexPrintf("COCO WARNING: ");
    va_start(args, message);
    n = vsprintf(buffer, message, args);
    va_end(args);
    mexPrintf(buffer);

    mexPrintf("\n");
  }
}

void coco_info(char const* message, ...) {
  va_list args;

  char buffer[500];
  int n;

  if (coco_log_level >= COCO_INFO) {
    mexPrintf("COCO INFO: ");
    va_start(args, message);
    n = vsprintf(buffer, message, args);
    va_end(args);
    mexPrintf(buffer);
    mexPrintf("\n");
  }
}

/**
 * A function similar to coco_info that prints only the given message without any prefix and without
 * adding a new line.
 */
void coco_info_partial(char const* message, ...) {
  va_list args;

  char buffer[500];
  int n;

  if (coco_log_level >= COCO_INFO) {
    va_start(args, message);
    n = vsprintf(buffer, message, args);
    va_end(args);
    mexPrintf(buffer);
  }
}

void coco_debug(char const* message, ...) {
  va_list args;

  char buffer[500];
  int n;

  if (coco_log_level >= COCO_DEBUG) {
    mexPrintf("COCO DEBUG: ");
    va_start(args, message);
    n = vsprintf(buffer, message, args);
    va_end(args);
    mexPrintf(buffer);
    mexPrintf("\n");
  }
}

void* coco_allocate_memory(const size_t size) {
  void* data;
  if (size == 0) {
    coco_error("coco_allocate_memory() called with 0 size.");
    return nullptr; /* never reached */
  }
  data = malloc(size);
  if (! data)
    coco_error("coco_allocate_memory() failed.");
  return data;
}

void coco_free_memory(void* data) { free(data); }
