/**
 * @file coco_utilities.cpp
 * @brief Definitions of miscellaneous functions used throughout the COCO framework.
 */

#include "coco_platform.h"

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>

#include <algorithm>

#include "coco.h"
#include "coco_internal.h"
#include "coco_string.cpp"

/***********************************************************************************************************/

/**
 * @brief Sets the constant chosen_precision to 1e-9.
 */
static const double chosen_precision = 1e-9;

/***********************************************************************************************************/

/**
 * @brief Initializes the logging level to COCO_INFO.
 */
static coco_log_level_type_e coco_log_level = COCO_INFO;

/**
 * @param log_level Denotes the level of information given to the user through the standard output and
 * error streams. Can take on the values:
 * - "error" (only error messages are output),
 * - "warning" (only error and warning messages are output),
 * - "info" (only error, warning and info messages are output) and
 * - "debug" (all messages are output).
 * - "" does not set a new value
 * The default value is info.
 *
 * @return The previous coco_log_level value as an immutable string.
 */
const char *coco_set_log_level(const char *log_level) {

  coco_log_level_type_e previous_log_level = coco_log_level;

  if (strcmp(log_level, "error") == 0)
    coco_log_level = COCO_ERROR;
  else if (strcmp(log_level, "warning") == 0)
    coco_log_level = COCO_WARNING;
  else if (strcmp(log_level, "info") == 0)
    coco_log_level = COCO_INFO;
  else if (strcmp(log_level, "debug") == 0)
    coco_log_level = COCO_DEBUG;
  else if (strcmp(log_level, "") == 0) {
    /* Do nothing */
  } else {
    coco_warning("coco_set_log_level(): unknown level %s", log_level);
  }

  if (previous_log_level == COCO_ERROR)
    return "error";
  else if (previous_log_level == COCO_WARNING)
    return "warning";
  else if (previous_log_level == COCO_INFO)
    return "info";
  else if (previous_log_level == COCO_DEBUG)
    return "debug";
  else {
    coco_error("coco_set_log_level(): unknown previous log level");
    return "";
  }
}

/***********************************************************************************************************/

/**
 * @name Methods regarding file, directory and path manipulations
 */
/**@{*/
/**
 * @brief Returns the given path (end of the recursion of the variadic coco_join_path()).
 */
static std::string coco_join_path(const std::string &path) { return path; }

/**
 * @brief Creates a platform-dependent path from the given strings.
 *
 * Usage examples:
 * - coco_join_path(base_path, folder1, folder2, folder3) returns base_path/folder1/folder2/folder3
 * - coco_join_path(base_path, folder1, file_name) returns base_path/folder1/file_name
 *
 * No separator is added after an empty base path.
 */
template <typename... Components>
static std::string coco_join_path(const std::string &path, const std::string &component,
                                  const Components &...components) {
  return coco_join_path(path.empty() ? component : path + coco_path_separator + component, components...);
}

/**
 * @brief Checks if the given directory exists.
 *
 * @note Should work cross-platform.
 *
 * @param path The given path.
 *
 * @return 1 if the path exists and corresponds to a directory and 0 otherwise.
 */
static int coco_directory_exists(const char *path) {
  int res;
#if defined(HAVE_GFA)
  DWORD dwAttrib = GetFileAttributesA(path);
  res = (dwAttrib != INVALID_FILE_ATTRIBUTES && (dwAttrib & FILE_ATTRIBUTE_DIRECTORY));
#elif defined(HAVE_STAT)
  struct stat buf;
  res = (!stat(path, &buf) && S_ISDIR(buf.st_mode));
#else
#error Ooops
#endif
  return res;
}

/**
 * @brief Checks if the given file exists.
 *
 * @note Should work cross-platform.
 *
 * @param path The given path.
 *
 * @return 1 if the path exists and corresponds to a file and 0 otherwise.
 */
static int coco_file_exists(const char *path) {
  int res;
#if defined(HAVE_GFA)
  DWORD dwAttrib = GetFileAttributesA(path);
  res = (dwAttrib != INVALID_FILE_ATTRIBUTES) && !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY);
#elif defined(HAVE_STAT)
  struct stat buf;
  res = (!stat(path, &buf) && !S_ISDIR(buf.st_mode));
#else
#error Ooops
#endif
  return res;
}

/**
 * @brief Calls the right mkdir() method (depending on the platform) with full privileges for the user.
 * If the created directory has not existed before, returns 0, otherwise returns 1. If the directory has
 * not been created, a coco_error is raised.
 *
 * @param path The directory path.
 *
 * @return 0 if the created directory has not existed before and 1 otherwise.
 */
static int coco_mkdir(const char *path) {
  int result = 0;

  /* Do not create the path if is of the form "C:" (two letters, of which the second is a colon)*/
  if ((strlen(path) == 2) && (path[1] == ':'))
    return 1;

#if _MSC_VER
  result = _mkdir(path);
#elif defined(__MINGW32__) || defined(__MINGW64__)
  result = mkdir(path);
#else
  result = mkdir(path, S_IRWXU);
#endif

  if (result == 0)
    return 0;
  else if (errno == EEXIST)
    return 1;
  else {
    coco_error("coco_mkdir(): unable to create %s, mkdir error: %s", path, strerror(errno));
    return 1; /* Never reached */
  }
}

/**
 * @brief Creates a directory (possibly having to create nested directories). If the last created directory
 * has not existed before, returns 0, otherwise returns 1.
 *
 * @param path The directory path.
 *
 * @return 0 if the created directory has not existed before and 1 otherwise.
 */
static int coco_create_directory(const std::string &path) {
  const char path_sep = coco_path_separator[0];
  std::string tmp = path;
  size_t i;

  /* Remove possible trailing (back)slash */
  if (!tmp.empty() && (tmp.back() == path_sep))
    tmp.pop_back();

  /* Iterate through nested directories (does nothing if directories are not nested) */
  for (i = 1; i < tmp.size(); i++) {
    if (tmp[i] == path_sep)
      coco_mkdir(tmp.substr(0, i).c_str());
  }

  /* Create the last nested or only directory */
  return coco_mkdir(tmp.c_str());
}

/**
 * @brief Returns a unique file name created from the given file_name.
 *
 * If path/file_name.ext does not yet exit, file_name is returned, otherwise a number is appended to it. If
 * path/file_name.ext already exists, path/filename-0001.ext will be tried. If this one exists as well,
 * path/filename-0002.ext will be tried, and so on. If path/filename-9999.ext exists as well, the function
 * raises an error. Every 1000 trials a warning is issued.
 */
static std::string coco_create_unique_filename(const std::string &path, const std::string &file_name,
                                               const char *ext) {
  int counter = 1;

  if (!coco_file_exists(coco_join_path(path, file_name + ext).c_str()))
    return file_name;

  while (counter < 9999) {
    const std::string new_file_name = coco_string_format("%s-%04d", file_name.c_str(), counter);

    if (!coco_file_exists(coco_join_path(path, new_file_name + ext).c_str()))
      return new_file_name;

    counter++;
    if (counter % 1000 == 0)
      coco_warning("coco_create_unique_filename(): trying to create a unique file name %s (%d trials)",
                   file_name.c_str(), counter);
  }

  coco_error("coco_create_unique_filename(): could not create a unique file name %s", file_name.c_str());
  return file_name; /* Never reached */
}

/**
 * @brief Creates a directory that has not existed before and returns its path.
 *
 * If the given path does not yet exit, it is used as is, otherwise a number is appended to it. If path
 * already exists, path-0001 will be tried. If this one exists as well, path-0002 will be tried, and so on.
 * If path-9999 exists as well, an error is raised. Every 1000 trials a warning is issued.
 */
static std::string coco_create_unique_directory(const std::string &path) {

  int counter = 1;

  if (coco_create_directory(path) == 0) {
    /* Directory created */
    return path;
  }

  while (counter < 9999) {
    const std::string new_path = coco_string_format("%s-%04d", path.c_str(), counter);

    if (coco_create_directory(new_path) == 0) {
      /* Directory created */
      return new_path;
    }
    counter++;
    if (counter % 1000 == 0)
      coco_warning("coco_create_unique_directory(): creating a unique directory %s (%d trials)", path.c_str(),
                   counter);
  }

  coco_error("coco_create_unique_directory(): unable to create unique directory %s", path.c_str());
  return path; /* Never reached */
}

/**
 * The method should work across different platforms/compilers.
 *
 * @path The path to the directory
 *
 * @return 0 on successful completion, and -1 on error.
 */
int coco_remove_directory(const char *path) {
#if _MSC_VER
  WIN32_FIND_DATA find_data_file;
  HANDLE find_handle = nullptr;
  int r = -1;
  int r2 = -1;

  /* Nothing to do if the folder does not exist */
  if ((find_handle = FindFirstFile(coco_string_format("%s\\*.*", path).c_str(), &find_data_file)) ==
      INVALID_HANDLE_VALUE) {
    return 0;
  }

  do {
    r = 0;

    /* Skip the names "." and ".." as we don't want to recurse on them */
    if (strcmp(find_data_file.cFileName, ".") != 0 && strcmp(find_data_file.cFileName, "..") != 0) {
      /* Build the new path using the argument path the file/folder name we just found */
      const std::string buf = coco_string_format("%s\\%s", path, find_data_file.cFileName);

      if (find_data_file.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
        /* Buf is a directory, recurse on it */
        r2 = coco_remove_directory(buf.c_str());
      } else {
        /* Buf is a file, delete it */
        /* Careful, DeleteFile returns 0 if it fails and nonzero otherwise! */
        r2 = -(DeleteFile(buf.c_str()) == 0);
      }
    }

    r = r2;

  } while (FindNextFile(find_handle, &find_data_file)); /* Find the next file */

  FindClose(find_handle);

  if (!r) {
    /* Path is an empty directory, delete it */
    /* Careful, RemoveDirectory returns 0 if it fails and nonzero otherwise! */
    r = -(RemoveDirectory(path) == 0);
  }

  return r;
#else
  DIR *d = opendir(path);
  int r = -1;
  int r2 = -1;

  /* Nothing to do if the folder does not exist */
  if (!coco_directory_exists(path))
    return 0;

  if (d) {
    struct dirent *p;

    r = 0;

    while (!r && (p = readdir(d))) {

      /* Skip the names "." and ".." as we don't want to recurse on them */
      if (!strcmp(p->d_name, ".") || !strcmp(p->d_name, "..")) {
        continue;
      }

      const std::string buf = coco_string_format("%s/%s", path, p->d_name);
      if (coco_directory_exists(buf.c_str())) {
        /* Buf is a directory, recurse on it */
        r2 = coco_remove_directory(buf.c_str());
      } else {
        /* Buf is a file, delete it */
        r2 = unlink(buf.c_str());
      }

      r = r2;
    }

    closedir(d);
  }

  if (!r) {
    /* Path is an empty directory, delete it */
    r = rmdir(path);
  }

  return r;
#endif
}

/**
 * The method should work across different platforms/compilers.
 *
 * @file_name The path to the file
 *
 * @return 0 on successful completion, and -1 on error.
 */
int coco_remove_file(const char *file_name) {
#if _MSC_VER
  int r = -1;
  /* Try to delete the file */
  /* Careful, DeleteFile returns 0 if it fails and nonzero otherwise! */
  r = -(DeleteFile(file_name) == 0);
  return r;
#else
  int r = -1;
  /* Try to delete the file */
  r = unlink(file_name);
  return r;
#endif
}
/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding memory allocations
 */
/**@{*/
double *coco_allocate_vector(const size_t number_of_elements) {
  const size_t block_size = number_of_elements * sizeof(double);
  return (double *)coco_allocate_memory(block_size);
}

/**
 * @brief A matrix of doubles stored row by row, whose rows may have different lengths.
 *
 * rows() returns an array of row pointers, as expected by the functions that work on double ** (such as
 * bbob2009_compute_rotation()). Block matrices only store the elements of the blocks, so their rows are
 * shorter than the number of columns. The elements are initialized to 0.
 *
 * The matrix can be moved but not copied, since copies would share the row pointers of the original.
 */
class coco_matrix_t {
public:
  coco_matrix_t() = default;

  /** @brief A dense matrix with the given number of rows and columns. */
  coco_matrix_t(const size_t number_of_rows, const size_t number_of_columns)
      : coco_matrix_t(std::vector<size_t>(number_of_rows, number_of_columns)) {}

  /** @brief A matrix whose row i has row_lengths[i] elements. */
  explicit coco_matrix_t(const std::vector<size_t> &row_lengths) : rows_(row_lengths.size()) {
    size_t i, offset = 0;
    for (i = 0; i < row_lengths.size(); ++i)
      offset += row_lengths[i];
    values_.resize(offset);
    for (i = 0, offset = 0; i < row_lengths.size(); ++i) {
      rows_[i] = values_.data() + offset;
      offset += row_lengths[i];
    }
  }

  coco_matrix_t(const coco_matrix_t &) = delete;
  coco_matrix_t &operator=(const coco_matrix_t &) = delete;
  coco_matrix_t(coco_matrix_t &&) = default;
  coco_matrix_t &operator=(coco_matrix_t &&) = default;

  double **rows() { return rows_.data(); }
  const double *const *rows() const { return rows_.data(); }
  double *operator[](const size_t i) { return rows_[i]; }
  const double *operator[](const size_t i) const { return rows_[i]; }
  size_t number_of_rows() const { return rows_.size(); }

private:
  std::vector<double> values_;
  std::vector<double *> rows_;
};
/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding string options
 */
/**@{*/

/**
 * @brief Returns redundant option keys (the ones present in given_option_keys but not in known_option_keys).
 */
static coco_option_keys_t coco_option_keys_get_redundant(const coco_option_keys_t &known_option_keys,
                                                         const coco_option_keys_t &given_option_keys) {
  coco_option_keys_t redundant_option_keys;

  for (const std::string &key : given_option_keys) {
    if (std::find(known_option_keys.begin(), known_option_keys.end(), key) == known_option_keys.end())
      redundant_option_keys.push_back(key);
  }
  return redundant_option_keys;
}

/**
 * @brief Returns the option keys from the given string of options containing keys and values separated by
 * colons.
 *
 * @note Relies heavily on the "key: value" format and might fail if the number of colons doesn't match the
 * number of keys. Values that are strings surrounded by quotation marks should work as long as they come
 * in pairs.
 */
static coco_option_keys_t coco_option_keys(const char *option_string) {

  size_t i;
  coco_option_keys_t option_keys;
  std::string cleaned_option_string;

  /* Check for empty string */
  if ((option_string == nullptr) || (strlen(option_string) == 0)) {
    return option_keys;
  }

  /* Construct the cleaned_option_string by replacing any string between two quotation marks with "STR" */
  const std::vector<std::string> parts = coco_string_split(option_string, '\"');
  for (i = 0; i < parts.size(); i++) {
    /* Even parts are outside of a pair of quotation marks */
    cleaned_option_string += (i % 2 == 0) ? parts[i] : "STR";
  }

  /* Split the options w.r.t ':' */
  const std::vector<std::string> keys = coco_string_split(cleaned_option_string, ':');

  /* Keys now contain something like this: "values_of_previous_key this_key" except for the first, which
   * contains only the key and the last, which contains only the previous values */
  for (i = 0; i < keys.size(); i++) {
    /* Remove any leading and trailing spaces */
    const std::string string_to_parse = coco_string_trim(keys[i]);

    /* Stop if this is the last substring (contains a value and no key) */
    if ((i > 0) && (i + 1 == keys.size()))
      break;

    /* Disregard everything before the last space (unless this is the first key) */
    const size_t space = string_to_parse.rfind(' ');
    if ((space == std::string::npos) || (i == 0))
      option_keys.push_back(string_to_parse);
    else
      option_keys.push_back(string_to_parse.substr(space + 1));
  }

  return option_keys;
}

/**
 * @brief Returns a string containing the info_string and all keys from option_keys (or an empty string if
 * there are no keys).
 *
 * Can be used to output information about the given option_keys.
 */
static std::string coco_option_keys_get_output_string(const coco_option_keys_t &option_keys,
                                                      const char *info_string) {
  std::string string;

  if (!option_keys.empty()) {
    string = info_string;
    for (const std::string &key : option_keys) {
      string += " " + key + "\n";
    }
  }
  return string;
}

/**
 * @brief Warns about the keys in options that are not among the known_keys.
 */
static void coco_option_keys_warn_redundant(const char *options, const coco_option_keys_t &known_option_keys,
                                            const char *redundant_info, const char *valid_info) {
  const coco_option_keys_t redundant_option_keys =
      coco_option_keys_get_redundant(known_option_keys, coco_option_keys(options));

  if (!redundant_option_keys.empty()) {
    /* Warn the user that some of given options are being ignored and output the valid options */
    coco_warning("%s%s", coco_option_keys_get_output_string(redundant_option_keys, redundant_info).c_str(),
                 coco_option_keys_get_output_string(known_option_keys, valid_info).c_str());
  }
}

/**
 * @brief Parses options in the form "name1: value1 name2: value2".
 *
 * Formatting requirements:
 * - name and value need to be separated by a colon (spaces are optional)
 * - value needs to be a single string (no spaces allowed)
 *
 * @return The number of successful assignments.
 */
static int coco_options_read(const char *options, const char *name, const char *format, void *pointer) {

  long i1, i2;

  if ((!options) || (strlen(options) == 0))
    return 0;

  i1 = coco_strfind(options, name);
  if (i1 < 0)
    return 0;
  i2 = i1 + coco_strfind(&options[i1], ":") + 1;

  /* Remove trailing whitespaces */
  while (isspace((unsigned char)options[i2]))
    i2++;

  if (i2 <= i1) {
    return 0;
  }

  return sscanf(&options[i2], format, pointer);
}

/**
 * @brief Reads an integer from options using the form "name1: value1 name2: value2".
 *
 * Formatting requirements:
 * - name and value need to be separated by a colon (spaces are optional)
 * - the value corresponding to the given name needs to be an integer
 *
 * @return The number of successful assignments.
 */
static int coco_options_read_int(const char *options, const char *name, int *pointer) {
  return coco_options_read(options, name, " %i", pointer);
}

/**
 * @brief Reads a size_t from options using the form "name1: value1 name2: value2".
 *
 * Formatting requirements:
 * - name and value need to be separated by a colon (spaces are optional)
 * - the value corresponding to the given name needs to be a size_t
 *
 * @return The number of successful assignments.
 */
static int coco_options_read_size_t(const char *options, const char *name, size_t *pointer) {
  return coco_options_read(options, name, "%lu", pointer);
}

/**
 * @brief Reads a double value from options using the form "name1: value1 name2: value2".
 *
 * Formatting requirements:
 * - name and value need to be separated by a colon (spaces are optional)
 * - the value corresponding to the given name needs to be a double
 *
 * @return The number of successful assignments.
 */
static int coco_options_read_double(const char *options, const char *name, double *pointer) {
  return coco_options_read(options, name, "%lf", pointer);
}

/**
 * @brief Returns the position of the value of the option with the given name in options (or -1 if the
 * option is not found).
 */
static long coco_options_find_value(const char *options, const char *name) {

  long i1, i2;

  if ((!options) || (strlen(options) == 0))
    return -1;

  i1 = coco_strfind(options, name);
  if (i1 < 0)
    return -1;
  i2 = i1 + coco_strfind(&options[i1], ":") + 1;

  if (i2 <= i1)
    return -1;
  return i2;
}

/**
 * @brief Reads a string from options using the form "name1: value1 name2: value2".
 *
 * Formatting requirements:
 * - name and value need to be separated by a colon (spaces are optional)
 * - the value corresponding to the given name needs to be a string - either a single word or multiple words
 * in double quotes
 *
 * @return The number of successful assignments (0 or 1). The value is left unchanged if nothing was read.
 */
static int coco_options_read_string(const char *options, const char *name, std::string &value) {

  const long i = coco_options_find_value(options, name);
  const char *start, *stop;

  if (i < 0)
    return 0;

  /* Remove leading white spaces */
  start = &options[i];
  while (isspace((unsigned char)*start))
    start++;

  if (*start == '\"') {
    /* The value starts with a quote: read everything between two quotes */
    start++;
    stop = strchr(start, '\"');
    if (stop == nullptr)
      stop = start + strlen(start);
  } else {
    stop = start;
    while ((*stop != '\0') && !isspace((unsigned char)*stop))
      stop++;
  }

  if (stop == start)
    return 0;
  value.assign(start, stop);
  return 1;
}

/**
 * @brief Reads (possibly delimited) values from options using the form "name1: value1,value2,value3 name2: value4",
 * i.e. reads all characters from the corresponding name up to the next alphabetic character or end of string,
 * ignoring white-space characters.
 *
 * Formatting requirements:
 * - names have to start with alphabetic characters
 * - values cannot include alphabetic characters
 * - name and value need to be separated by a colon (spaces are optional)
 *
 * @return The number of characters read.
 */
static size_t coco_options_read_values(const char *options, const char *name, std::string &value) {

  const long i = coco_options_find_value(options, name);
  const char *ptr;

  value.clear();
  if (i < 0)
    return 0;

  for (ptr = &options[i]; (*ptr != '\0') && !isalpha((unsigned char)*ptr); ptr++) {
    if (!isspace((unsigned char)*ptr))
      value += *ptr;
  }
  return value.size();
}
/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods implementing functions on double values not contained in C89 standard
 */
/**@{*/

/**
 * @brief  Returns 1 if |a - b| < precision and 0 otherwise.
 */
static int coco_double_almost_equal(const double a, const double b, const double precision) {
  return (fabs(a - b) < precision);
}

/**
 * @brief Rounds the given double to the nearest integer.
 */
static double coco_double_round(const double number) { return floor(number + 0.5); }

/**
 * @brief Rounds the given double up to the nearest double value with the given precision.
 *
 * @note The implementation is (probably unnecessarily) complex, but this was the only way to make
 * sure it works also for edge cases due to float precision issue.
 */
static double coco_double_round_up_with_precision(const double number, const double precision) {
  double rounded_up, rounded;
  double min_precision = 1e-12;
  assert(precision > min_precision);
  rounded_up = ceil(number / precision) * precision;
  rounded = coco_double_round(number / precision) * precision;
  if (coco_double_almost_equal(rounded, rounded_up, precision))
    return rounded_up;
  else {
    if (coco_double_almost_equal(number - rounded, 0, min_precision))
      return rounded;
    else
      return rounded_up;
  }
}

/**
 * @brief Returns the maximum of a and b.
 */
static double coco_double_max(const double a, const double b) {
  if (a >= b) {
    return a;
  } else {
    return b;
  }
}

/**
 * @brief Returns the minimum of a and b.
 */
static double coco_double_min(const double a, const double b) {
  if (a <= b) {
    return a;
  } else {
    return b;
  }
}

/**
 * @brief Performs a "safer" double to size_t conversion.
 *
 * TODO: This method could (should?) check for overflow when casting (similarly as is done in
 * coco_double_to_int()).
 */
static size_t coco_double_to_size_t(const double number) { return (size_t)coco_double_round(number); }

/**
 * @brief Rounds the given double to the nearest integer (returns the number in int type)
 */
static int coco_double_to_int(const double number) {
  if (number > (double)INT_MAX) {
    coco_error("coco_double_to_int(): Cannot cast %f to the nearest integer, max %d allowed", number, INT_MAX);
    return -1; /* Never reached */
  } else if (number < (double)INT_MIN) {
    coco_error("coco_double_to_int(): Cannot cast %f to the nearest integer, min %d allowed", number, INT_MIN);
    return -1; /* Never reached */
  } else {
    return (int)(number + 0.5);
  }
}

/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods handling NAN and INFINITY
 */
/**@{*/

/**
 * @brief Returns 1 if x is NAN and 0 otherwise.
 */
static int coco_is_nan(const double x) {
  return (isnan(x) || (x != x) || !(x == x) ||
          ((x >= NAN / (1 + chosen_precision)) && (x <= NAN * (1 + chosen_precision))));
}

/**
 * @brief Returns 1 if the input vector of dimension dim contains any NAN values and 0 otherwise.
 */
static int coco_vector_contains_nan(const double *x, const size_t dim) {
  size_t i;
  for (i = 0; i < dim; i++) {
    if (COCO_UNLIKELY(coco_is_nan(x[i])))
      return 1;
  }
  return 0;
}

/**
 * @brief Sets all dim values of y to NAN.
 */
static void coco_vector_set_to_nan(double *y, const size_t dim) {
  size_t i;
  for (i = 0; i < dim; i++) {
    y[i] = NAN;
  }
}

/**
 * @brief Returns 1 if x is INFINITY and 0 otherwise.
 */
static int coco_is_inf(const double x) {
  if (coco_is_nan(x))
    return 0;
  return (isinf(x) || (x <= -INFINITY) || (x >= INFINITY));
}

/**
 * @brief Sets all dim values of y to Inf.
 */
static void coco_vector_set_to_inf(double *y, const size_t dim) {
  size_t i;
  for (i = 0; i < dim; i++) {
    y[i] = INFINITY;
  }
}

/**
 * @brief Returns 1 if the input vector of dimension dim contains no NaN of inf values, and 0 otherwise.
 */
static int coco_vector_isfinite(const double *x, const size_t dim) {
  size_t i;
  for (i = 0; i < dim; i++) {
    if (COCO_UNLIKELY(coco_is_nan(x[i])) || COCO_UNLIKELY(coco_is_inf(x[i])))
      return 0;
  }
  return 1;
}

/**
 * @brief Returns 1 if the point x is feasible, and 0 otherwise.
 *
 * Allows constraint_values == NULL, otherwise constraint_values
 * must be a valid double* pointer and contains the g-values of x
 * on "return".
 *
 * Any point x containing NaN or inf values is considered infeasible.
 *
 * This function is (and should be) used internally only, and does not
 * increase the counter of constraint function evaluations.
 *
 * @param problem The given COCO problem.
 * @param x Decision vector.
 * @param constraint_values Vector of contraints values resulting from evaluation.
 */
static int coco_is_feasible(coco_problem_t *problem, const double *x, double *constraint_values) {

  size_t i;
  double *cons_values = constraint_values;
  std::vector<double> buffer;
  int ret_val = 1;

  /* Return 0 if the decision vector contains any INFINITY or NaN values */
  if (!coco_vector_isfinite(x, coco_problem_get_dimension(problem)))
    return 0;

  if (coco_problem_get_number_of_constraints(problem) <= 0)
    return 1;

  assert(problem != nullptr);

  if (constraint_values == nullptr) {
    buffer.resize(problem->number_of_constraints);
    cons_values = buffer.data();
  }

  problem->evaluate_constraint(x, cons_values, 0);

  for (i = 0; i < coco_problem_get_number_of_constraints(problem); ++i) {
    if (cons_values[i] > 0.0) {
      ret_val = 0;
      break;
    }
  }

  return ret_val;
}

/**@}*/

/***********************************************************************************************************/

/**
 * @name Miscellaneous methods
 */
/**@{*/

/**
 * @brief Returns the current time as a string.
 */
static std::string coco_current_time_get_string(void) {
  time_t timer;
  char time_string[30];
  struct tm *tm_info;
  time(&timer);
  tm_info = localtime(&timer);
  assert(tm_info != nullptr);
  strftime(time_string, sizeof(time_string), "%d.%m.%y %H:%M:%S", tm_info);
  return time_string;
}

/**
 * @brief multiply each componenent by nom/denom or by nom if denom == 0.
 *
 * return used scaling factor, usually nom/denom.
 *
 * Example: coco_vector_scale(x, dimension, 1, coco_vector_norm(x, dimension));
 */
static double coco_vector_scale(double *x, size_t dimension, double nom, double denom) {

  size_t i;

  assert(x);

  if (denom != 0)
    nom /= denom;

  for (i = 0; i < dimension; ++i)
    x[i] *= nom;
  return nom;
}

/**
 * @brief return norm of vector x.
 *
 */
static double coco_vector_norm(const double *x, size_t dimension) {

  size_t i;
  double ssum = 0.0;

  assert(x);

  for (i = 0; i < dimension; ++i)
    ssum += x[i] * x[i];

  return sqrt(ssum);
}

/**
 * @brief return scalar product between vectors x and y.
 *
 */
static double coco_vector_scalar_product(const double *x, const double *y, size_t dimension) {

  size_t i;
  double ssum = 0.0;

  assert(x);
  assert(y);

  for (i = 0; i < dimension; ++i)
    ssum += x[i] * y[i];

  return ssum;
}

/**
 * @brief Checks if a given matrix M is orthogonal by (partially) computing M * M^T.
 * If M is a square matrix and M * M^T is close enough to the identity matrix
 * (up to a chosen precision), the function returns 1. Otherwise, it returns 0.
 * The matrix M must be represented as an array of doubles.
 */
static int coco_is_orthogonal(const double *M, const size_t nb_rows, const size_t nb_columns) {

  size_t i, j, z;
  double sum;

  if (nb_rows != nb_columns)
    return 0;

  for (i = 0; i < nb_rows; ++i) {
    for (j = 0; j < nb_rows; ++j) {
      /* Compute the dot product of the ith row of M
       * and the jth column of M^T (i.e. jth row of M)
       */
      sum = 0.0;
      for (z = 0; z < nb_rows; ++z) {
        sum += M[i * nb_rows + z] * M[j * nb_rows + z];
      }

      /* Check if the dot product is 1 (resp. 0) when the row and the column
       * indices are the same (resp. different)
       */
      if (((i == j) && !coco_double_almost_equal(sum, 1, chosen_precision)) ||
          ((i != j) && !coco_double_almost_equal(sum, 0, chosen_precision)))
        return 0;
    }
  }
  return 1;
}

/**
 * @brief Returns 1 if the input vector x is (close to) zero and 0 otherwise.
 */
static int coco_vector_is_zero(const double *x, const size_t dim) {
  size_t i = 0;
  int is_zero = 1;

  if (coco_vector_contains_nan(x, dim))
    return 0;

  while (i < dim && is_zero) {
    is_zero = coco_double_almost_equal(x[i], 0, chosen_precision);
    i++;
  }

  return is_zero;
}
/**@}*/

/***********************************************************************************************************/
