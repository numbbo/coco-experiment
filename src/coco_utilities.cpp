/**
 * @file coco_utilities.cpp
 * @brief Definitions of miscellaneous functions used throughout the COCO framework.
 */

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>

#if __cplusplus < 202002L && (!defined(_MSVC_LANG) || _MSVC_LANG < 202002L)
#error "COCO needs C++20 (for example -std=c++20 or /std:c++20)"
#endif

#include <algorithm>
#include <filesystem>
#include <system_error>

/* Silence the Visual Studio warnings (C4996) about fopen(), sscanf() and similar functions */
#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

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
const char *coco_set_log_level(const char *log_level_string) {

  const std::string_view log_level = coco_c_string_view(log_level_string);
  coco_log_level_type_e previous_log_level = coco_log_level;

  if (log_level == "error")
    coco_log_level = COCO_ERROR;
  else if (log_level == "warning")
    coco_log_level = COCO_WARNING;
  else if (log_level == "info")
    coco_log_level = COCO_INFO;
  else if (log_level == "debug")
    coco_log_level = COCO_DEBUG;
  else if (log_level.empty()) {
    /* Do nothing */
  } else {
    coco_warning("coco_set_log_level(): unknown level %s", log_level_string);
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
 * @brief Creates a path from the given strings, using the platform's preferred separator.
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
  return coco_join_path(path.empty() ? component : (std::filesystem::path(path) / component).string(), components...);
}

/**
 * @brief Returns whether the path exists and is a directory.
 */
static bool coco_directory_exists(const std::string &path) {
  std::error_code error;
  return std::filesystem::is_directory(path, error);
}

/**
 * @brief Returns whether the path exists and is not a directory.
 */
static bool coco_file_exists(const std::string &path) {
  std::error_code error;
  const std::filesystem::file_status status = std::filesystem::status(path, error);
  return std::filesystem::exists(status) && !std::filesystem::is_directory(status);
}

/**
 * @brief Creates a directory, including any missing parent directories.
 *
 * Raises a coco_error if the directory cannot be created, for example because a file of the same name exists.
 *
 * @return Whether the directory was created (false if it existed before).
 */
static bool coco_create_directory(const std::string &path) {
  std::string directory = path;
  std::error_code error;

  /* Remove trailing separators, some implementations report a path like "a/b/" as not created */
  while (directory.size() > 1 &&
         (directory.back() == '/' || directory.back() == static_cast<char>(std::filesystem::path::preferred_separator)))
    directory.pop_back();

  const bool created = std::filesystem::create_directories(directory, error);
  if (error)
    coco_error("coco_create_directory(): unable to create %s: %s", path.c_str(), error.message().c_str());
  return created;
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

  if (coco_create_directory(path)) {
    /* Directory created */
    return path;
  }

  while (counter < 9999) {
    const std::string new_path = coco_string_format("%s-%04d", path.c_str(), counter);

    if (coco_create_directory(new_path)) {
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
 * @brief Removes the directory with all its contents. Does nothing if the directory does not exist.
 *
 * @param path The path to the directory.
 *
 * @return 0 on successful completion, and -1 on error.
 */
int coco_remove_directory(const char *path) {
  std::error_code error;

  if (!coco_directory_exists(path))
    return 0;
  std::filesystem::remove_all(path, error);
  return error ? -1 : 0;
}
/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding memory allocations
 */
/**@{*/
double *coco_allocate_vector(const size_t number_of_elements) {
  const size_t block_size = number_of_elements * sizeof(double);
  return static_cast<double *>(coco_allocate_memory(block_size));
}

/**
 * @brief A dense matrix of doubles, stored row by row (row-major) in a single vector.
 *
 * The elements are initialized to 0.
 */
class coco_matrix_t {
public:
  coco_matrix_t() = default;

  /** @brief A matrix with the given number of rows and columns. */
  coco_matrix_t(const size_t number_of_rows, const size_t number_of_columns)
      : number_of_rows_(number_of_rows), number_of_columns_(number_of_columns),
        values_(number_of_rows * number_of_columns) {}

  size_t number_of_rows() const { return number_of_rows_; }
  size_t number_of_columns() const { return number_of_columns_; }

  /** @brief The element in row i and column j. */
  double &operator()(const size_t i, const size_t j) { return values_[i * number_of_columns_ + j]; }
  double operator()(const size_t i, const size_t j) const { return values_[i * number_of_columns_ + j]; }

  /** @brief The elements of row i. */
  std::span<double> row(const size_t i) { return {values_.data() + i * number_of_columns_, number_of_columns_}; }
  std::span<const double> row(const size_t i) const {
    return {values_.data() + i * number_of_columns_, number_of_columns_};
  }

  /** @brief Returns whether the matrix is square and M M^T is the identity up to chosen_precision. */
  bool is_orthogonal() const;

  /**
   * @brief Solves M x = y by Gaussian elimination with partial pivoting. Returns false if the matrix is not
   * square or singular.
   */
  bool solve(std::span<const double> y, std::span<double> x) const;

  /** @brief Computes y = M x + b; x and y must not overlap (b and y may be the same). */
  void multiply_add(std::span<const double> x, std::span<const double> b, std::span<double> y) const {
    assert(x.size() == number_of_columns_ && b.size() == number_of_rows_ && y.size() == number_of_rows_);
    for (size_t i = 0; i < number_of_rows_; ++i) {
      const std::span<const double> current_row = row(i);
      y[i] = b[i];
      for (size_t j = 0; j < number_of_columns_; ++j)
        y[i] += current_row[j] * x[j];
    }
  }

  /** @brief Computes y = M x; x and y must not overlap. */
  void multiply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == number_of_columns_ && y.size() == number_of_rows_);
    for (size_t i = 0; i < number_of_rows_; ++i) {
      const std::span<const double> current_row = row(i);
      y[i] = 0;
      for (size_t j = 0; j < number_of_columns_; ++j)
        y[i] += current_row[j] * x[j];
    }
  }

  /** @brief Computes y = M^T x; x and y must not overlap. */
  void multiply_transposed(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == number_of_rows_ && y.size() == number_of_columns_);
    for (size_t j = 0; j < number_of_columns_; ++j)
      y[j] = 0;
    for (size_t i = 0; i < number_of_rows_; ++i) {
      const std::span<const double> current_row = row(i);
      for (size_t j = 0; j < number_of_columns_; ++j)
        y[j] += current_row[j] * x[i];
    }
  }

private:
  size_t number_of_rows_ = 0;
  size_t number_of_columns_ = 0;
  std::vector<double> values_;
};

/**
 * @brief Returns A diag(d) B, with every element computed as the sum over k of (A(i, k) * d[k]) * B(k, j) in the
 * order of k.
 */
static coco_matrix_t coco_matrix_product_with_diagonal(const coco_matrix_t &A, const std::vector<double> &d,
                                                       const coco_matrix_t &B) {
  coco_matrix_t result(A.number_of_rows(), B.number_of_columns());
  assert(A.number_of_columns() == d.size() && B.number_of_rows() == d.size());
  for (size_t i = 0; i < result.number_of_rows(); ++i) {
    for (size_t j = 0; j < result.number_of_columns(); ++j) {
      for (size_t k = 0; k < d.size(); ++k)
        result(i, j) += A(i, k) * d[k] * B(k, j);
    }
  }
  return result;
}

/**
 * @brief A square block-diagonal matrix of which only the blocks on the diagonal are stored.
 *
 * Row i holds the block_size(i) elements of its block, which lie in the columns first_column(i), ...,
 * first_column(i) + block_size(i) - 1; all other elements are 0. The elements are initialized to 0.
 */
class coco_block_diagonal_matrix_t {
public:
  coco_block_diagonal_matrix_t() = default;

  /** @brief A matrix with number_of_blocks blocks of the given sizes. */
  explicit coco_block_diagonal_matrix_t(std::span<const size_t> block_sizes) {
    size_t first_column = 0, number_of_values = 0;
    for (const size_t block_size : block_sizes) {
      for (size_t r = 0; r < block_size; ++r) {
        block_size_.push_back(block_size);
        first_column_.push_back(first_column);
        row_start_.push_back(number_of_values);
        number_of_values += block_size;
      }
      first_column += block_size;
    }
    values_.resize(number_of_values);
  }

  /** @brief The number of rows (and columns). */
  size_t dimension() const { return block_size_.size(); }
  /** @brief The size of the block that row i belongs to. */
  size_t block_size(const size_t i) const { return block_size_[i]; }
  /** @brief The first column of the block that row i belongs to. */
  size_t first_column(const size_t i) const { return first_column_[i]; }

  /** @brief The element in row i and column first_column(i) + k. */
  double &in_block(const size_t i, const size_t k) { return values_[row_start_[i] + k]; }
  double in_block(const size_t i, const size_t k) const { return values_[row_start_[i] + k]; }

  /** @brief The element in row i and column j (0 outside the blocks). */
  double operator()(const size_t i, const size_t j) const {
    return (j < first_column_[i] || j >= first_column_[i] + block_size_[i]) ? 0 : in_block(i, j - first_column_[i]);
  }

  /** @brief Computes y = B x; x and y must not overlap. */
  void multiply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == dimension() && y.size() == dimension());
    for (size_t i = 0; i < dimension(); ++i) {
      y[i] = 0;
      for (size_t j = first_column_[i]; j < first_column_[i] + block_size_[i]; ++j)
        y[i] += in_block(i, j - first_column_[i]) * x[j];
    }
  }

  /** @brief Computes y = B^T x; x and y must not overlap. */
  void multiply_transposed(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == dimension() && y.size() == dimension());
    for (size_t j = 0; j < dimension(); ++j)
      y[j] = 0;
    for (size_t i = 0; i < dimension(); ++i) {
      for (size_t j = first_column_[i]; j < first_column_[i] + block_size_[i]; ++j)
        y[j] += in_block(i, j - first_column_[i]) * x[i];
    }
  }

private:
  std::vector<double> values_;
  std::vector<size_t> block_size_;
  std::vector<size_t> first_column_;
  std::vector<size_t> row_start_;
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
static coco_option_keys_t coco_option_keys(std::string_view option_string) {

  coco_option_keys_t option_keys;
  std::string cleaned_option_string;

  /* Check for empty string */
  if (option_string.empty()) {
    return option_keys;
  }

  /* Construct the cleaned_option_string by replacing any string between two quotation marks with "STR" */
  const std::vector<std::string> parts = coco_string_split(option_string, '\"');
  for (size_t i = 0; i < parts.size(); i++) {
    /* Even parts are outside of a pair of quotation marks */
    cleaned_option_string += (i % 2 == 0) ? parts[i] : "STR";
  }

  /* Split the options w.r.t ':' */
  const std::vector<std::string> keys = coco_string_split(cleaned_option_string, ':');

  /* Keys now contain something like this: "values_of_previous_key this_key" except for the first, which
   * contains only the key and the last, which contains only the previous values */
  for (size_t i = 0; i < keys.size(); i++) {
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
static void coco_option_keys_warn_redundant(std::string_view options, const coco_option_keys_t &known_option_keys,
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
 * @brief Returns the position of the value of the option with the given name in options (or npos if the
 * option is not found).
 */
static size_t coco_options_find_value(std::string_view options, std::string_view name) {

  const size_t name_position = options.find(name);
  if (name_position == std::string_view::npos)
    return std::string_view::npos;
  const size_t colon_position = options.find(':', name_position);
  if (colon_position == std::string_view::npos)
    return std::string_view::npos;
  return colon_position + 1;
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
static int coco_options_read(std::string_view options, std::string_view name, const char *format, void *pointer) {

  size_t start = coco_options_find_value(options, name);
  if (start == std::string_view::npos)
    return 0;

  /* Skip leading white space and take the value up to the next white space */
  while ((start < options.size()) && isspace(static_cast<unsigned char>(options[start])))
    start++;
  size_t stop = start;
  while ((stop < options.size()) && !isspace(static_cast<unsigned char>(options[stop])))
    stop++;

  const std::string value(options.substr(start, stop - start));
  return sscanf(value.c_str(), format, pointer);
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
static int coco_options_read_int(std::string_view options, std::string_view name, int *pointer) {
  return coco_options_read(options, name, " %i", pointer);
}

/**
 * @brief Reads a flag (an integer, nonzero meaning true) from options using the form "name1: value1 name2: value2".
 *
 * @return The number of successful assignments. The flag is left unchanged if nothing was read.
 */
static int coco_options_read_flag(std::string_view options, std::string_view name, bool *pointer) {
  int value = 0;
  const int result = coco_options_read_int(options, name, &value);
  if (result == 1)
    *pointer = (value != 0);
  return result;
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
static int coco_options_read_size_t(std::string_view options, std::string_view name, size_t *pointer) {
  return coco_options_read(options, name, "%zu", pointer);
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
static int coco_options_read_double(std::string_view options, std::string_view name, double *pointer) {
  return coco_options_read(options, name, "%lf", pointer);
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
static int coco_options_read_string(std::string_view options, std::string_view name, std::string &value) {

  size_t start = coco_options_find_value(options, name);
  size_t stop;

  if (start == std::string_view::npos)
    return 0;

  /* Remove leading white spaces */
  while ((start < options.size()) && isspace(static_cast<unsigned char>(options[start])))
    start++;

  if ((start < options.size()) && (options[start] == '\"')) {
    /* The value starts with a quote: read everything between two quotes */
    start++;
    stop = options.find('\"', start);
    if (stop == std::string_view::npos)
      stop = options.size();
  } else {
    stop = start;
    while ((stop < options.size()) && !isspace(static_cast<unsigned char>(options[stop])))
      stop++;
  }

  if (stop == start)
    return 0;
  value.assign(options.substr(start, stop - start));
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
static size_t coco_options_read_values(std::string_view options, std::string_view name, std::string &value) {

  value.clear();
  const size_t start = coco_options_find_value(options, name);
  if (start == std::string_view::npos)
    return 0;

  for (size_t i = start; (i < options.size()) && !isalpha(static_cast<unsigned char>(options[i])); i++) {
    if (!isspace(static_cast<unsigned char>(options[i])))
      value += options[i];
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
 * @brief  Returns whether |a - b| < precision.
 */
static bool coco_double_almost_equal(const double a, const double b, const double precision) {
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
static size_t coco_double_to_size_t(const double number) { return static_cast<size_t>(coco_double_round(number)); }

/**
 * @brief Rounds the given double to the nearest integer (returns the number in int type)
 */
static int coco_double_to_int(const double number) {
  if (number > static_cast<double>(INT_MAX)) {
    coco_error("coco_double_to_int(): Cannot cast %f to the nearest integer, max %d allowed", number, INT_MAX);
    return -1; /* Never reached */
  } else if (number < static_cast<double> INT_MIN) {
    coco_error("coco_double_to_int(): Cannot cast %f to the nearest integer, min %d allowed", number, INT_MIN);
    return -1; /* Never reached */
  } else {
    return static_cast<int>(number + 0.5);
  }
}

/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods handling NAN and INFINITY
 */
/**@{*/

/**
 * @brief Returns whether x is NAN.
 */
static bool coco_is_nan(const double x) {
  return (isnan(x) || (x != x) || !(x == x) ||
          ((x >= NAN / (1 + chosen_precision)) && (x <= NAN * (1 + chosen_precision))));
}

/**
 * @brief Returns whether the vector x contains any NAN values.
 */
static bool coco_vector_contains_nan(std::span<const double> x) {
  for (const double xi : x) {
    if (COCO_UNLIKELY(coco_is_nan(xi)))
      return true;
  }
  return false;
}

/**
 * @brief Sets all values of y to NAN.
 */
static void coco_vector_set_to_nan(std::span<double> y) {
  for (double &yi : y) {
    yi = NAN;
  }
}

/**
 * @brief Returns whether x is INFINITY.
 */
static bool coco_is_inf(const double x) {
  if (coco_is_nan(x))
    return false;
  return (isinf(x) || (x <= -INFINITY) || (x >= INFINITY));
}

/**
 * @brief Sets all values of y to Inf.
 */
static void coco_vector_set_to_inf(std::span<double> y) {
  for (double &yi : y) {
    yi = INFINITY;
  }
}

/**
 * @brief Returns whether the vector x contains no NaN or inf values.
 */
static bool coco_vector_isfinite(std::span<const double> x) {
  for (const double xi : x) {
    if (COCO_UNLIKELY(coco_is_nan(xi)) || COCO_UNLIKELY(coco_is_inf(xi)))
      return false;
  }
  return true;
}

/**
 * @brief Returns whether the point x is feasible.
 *
 * Allows an empty constraint_values, otherwise constraint_values
 * must have one element per constraint and contains the g-values of x
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
bool coco_problem_s::is_feasible(std::span<const double> x, std::span<double> constraint_values) {

  std::span<double> cons_values = constraint_values;
  std::vector<double> buffer;
  int ret_val = 1;

  /* Return 0 if the decision vector contains any INFINITY or NaN values */
  if (!coco_vector_isfinite(x))
    return false;

  if (coco_problem_get_number_of_constraints(this) <= 0)
    return true;

  if (constraint_values.empty()) {
    buffer.resize(number_of_constraints);
    cons_values = buffer;
  }

  evaluate_constraint(x, cons_values, 0);

  for (size_t i = 0; i < coco_problem_get_number_of_constraints(this); ++i) {
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
static std::string coco_current_time_get_string() {
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
 * Example: coco_vector_scale(x, 1, coco_vector_norm(x));
 */
static double coco_vector_scale(std::span<double> x, double nom, double denom) {
  if (denom != 0)
    nom /= denom;

  for (double &xi : x)
    xi *= nom;
  return nom;
}

/**
 * @brief return norm of vector x.
 *
 */
static double coco_vector_norm(std::span<const double> x) {

  double ssum = 0.0;

  for (const double xi : x)
    ssum += xi * xi;

  return sqrt(ssum);
}

/**
 * @brief return scalar product between vectors x and y.
 *
 */
static double coco_vector_scalar_product(std::span<const double> x, std::span<const double> y) {

  double ssum = 0.0;

  assert(x.size() == y.size());
  for (size_t i = 0; i < x.size(); ++i)
    ssum += x[i] * y[i];

  return ssum;
}

/**
 * @brief Checks whether the matrix is orthogonal by computing M * M^T. Returns true if M is square and M * M^T is
 * close enough to the identity matrix (up to chosen_precision).
 */
bool coco_matrix_t::is_orthogonal() const {
  double sum;

  if (number_of_rows_ != number_of_columns_)
    return false;

  for (size_t i = 0; i < number_of_rows_; ++i) {
    for (size_t j = 0; j < number_of_rows_; ++j) {
      /* Compute the dot product of the ith row of M and the jth column of M^T (i.e. jth row of M) */
      sum = 0.0;
      for (size_t z = 0; z < number_of_columns_; ++z) {
        sum += (*this)(i, z) * (*this)(j, z);
      }

      /* Check if the dot product is 1 (resp. 0) when the row and the column
       * indices are the same (resp. different)
       */
      if (((i == j) && !coco_double_almost_equal(sum, 1, chosen_precision)) ||
          ((i != j) && !coco_double_almost_equal(sum, 0, chosen_precision)))
        return false;
    }
  }
  return true;
}

bool coco_matrix_t::solve(std::span<const double> y, std::span<double> x) const {
  const size_t n = number_of_rows_;
  if (n != number_of_columns_)
    return false;

  assert(y.size() == n && x.size() == n);
  coco_matrix_t A = *this;
  std::vector<double> b(y.begin(), y.end());
  for (size_t k = 0; k < n; ++k) {
    size_t pivot = k;
    for (size_t i = k + 1; i < n; ++i) {
      if (fabs(A(i, k)) > fabs(A(pivot, k)))
        pivot = i;
    }
    if (A(pivot, k) == 0.0)
      return false;
    if (pivot != k) {
      for (size_t j = k; j < n; ++j)
        std::swap(A(k, j), A(pivot, j));
      std::swap(b[k], b[pivot]);
    }
    for (size_t i = k + 1; i < n; ++i) {
      const double factor = A(i, k) / A(k, k);
      for (size_t j = k; j < n; ++j)
        A(i, j) -= factor * A(k, j);
      b[i] -= factor * b[k];
    }
  }
  for (size_t i = n; i-- > 0;) {
    double sum = b[i];
    for (size_t j = i + 1; j < n; ++j)
      sum -= A(i, j) * x[j];
    x[i] = sum / A(i, i);
  }
  return true;
}

/**
 * @brief Returns whether the input vector x is (close to) zero.
 */
static bool coco_vector_is_zero(std::span<const double> x) {
  size_t i = 0;
  bool is_zero = true;

  if (coco_vector_contains_nan(x))
    return false;

  while (i < x.size() && is_zero) {
    is_zero = coco_double_almost_equal(x[i], 0, chosen_precision);
    i++;
  }

  return is_zero;
}
/**@}*/

/***********************************************************************************************************/
