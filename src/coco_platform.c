#include "coco_platform.h"

#if defined(_WIN32) || defined(_WIN64) || defined(__MINGW64__) || defined(__CYGWIN__)
const char *coco_path_separator = "\\";
#elif defined(__gnu_linux__) || \
      defined(__linux__) || \
      defined(__EMSCRIPTEN__) || \
      defined(__APPLE__) || \
      defined(__FreeBSD__) || \
      (defined(__sun) || defined(sun)) && (defined(__SVR4) || defined(__svr4__))
const char *coco_path_separator = "/";
#else
#error Unknown platform
#endif
