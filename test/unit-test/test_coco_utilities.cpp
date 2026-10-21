#include "minunit.h"

#include "coco.cpp"

static std::string convert_to_string_with_newlines(const std::vector<std::string> &array) {

  std::string result;

  for (const std::string &item : array) {
    result += item + "\n";
  }
  return result;
}

/**
 * Tests the function coco_set_log_level.
 */
MU_TEST(test_coco_set_log_level) {

  /* Remember the log level before all these tests */
  const std::string previous_log_level = coco_set_log_level("info");

  /* Check whether the level is set to COCO_INFO */
  mu_check(std::string(coco_set_log_level("")) == "info");

  /* Check whether the method works */
  mu_check(std::string(coco_set_log_level("error")) == "info");
  mu_check(std::string(coco_set_log_level("")) == "error");

  mu_check(std::string(coco_set_log_level("warning")) == "error");
  mu_check(std::string(coco_set_log_level("")) == "warning");

  mu_check(std::string(coco_set_log_level("debug")) == "warning");
  mu_check(std::string(coco_set_log_level("")) == "debug");

  mu_check(std::string(coco_set_log_level("info")) == "debug");
  mu_check(std::string(coco_set_log_level("")) == "info");

  /* An invalid argument shouldn't change the current value */
  coco_set_log_level("error");
  mu_check(std::string(coco_set_log_level("bla")) == "error");
  mu_check(std::string(coco_set_log_level("")) == "error");

  /* Restore the log level */
  coco_set_log_level(previous_log_level.c_str());
}

/**
 * Returns the prefix followed by the current local time, so that test directories do not collide.
 */
static std::string time_stamped(const std::string &prefix) {
  const time_t date_time = time(nullptr);
  char date_time_string[32];

  if (strftime(date_time_string, sizeof(date_time_string), "_%Y_%m_%d_%H_%M_%S", localtime(&date_time)) == 0)
    coco_error("time_stamped(): cannot decode local time string");
  return prefix + date_time_string;
}

/*
 * Tests creating and removing directory.
 */
MU_TEST(test_coco_create_remove_directory) {

  const std::string path = time_stamped("temp");

  /* At the beginning the path should not exist. */
  mu_check(!coco_directory_exists(path));

  coco_create_directory(path);
  mu_check(coco_directory_exists(path));

  /* Calling it again to check the handling if the path does exist. */
  coco_create_directory(path);

  coco_remove_directory(path.c_str());
  mu_check(!coco_directory_exists(path));

  /* Calling it again to check the handling if the path does not exist. */
  coco_remove_directory(path.c_str());
}

/*
 * Tests nested directories, unique names, file detection and removing a directory with contents.
 */
MU_TEST(test_coco_directories_and_files) {

  const std::string base = time_stamped("temp_files");
  const std::string nested = coco_join_path(base, "a", "b");

  mu_check(!coco_directory_exists(base));
  mu_check(coco_create_directory(nested)); /* Creates base, a and b */
  mu_check(coco_directory_exists(coco_join_path(base, "a")));
  mu_check(!coco_create_directory(nested)); /* Already exists */
  mu_check(!coco_create_directory(nested + "/"));

  /* A unique directory gets a numbered suffix when the name is taken */
  const std::string unique = coco_create_unique_directory(nested);
  mu_check(unique == nested + "-0001");
  mu_check(coco_directory_exists(unique));

  /* Files are found, directories are not files */
  const std::string file_name = coco_join_path(nested, "data.txt");
  mu_check(!coco_file_exists(file_name));
  FILE *file = fopen(file_name.c_str(), "w");
  mu_check(file != NULL);
  fprintf(file, "1\n");
  fclose(file);
  mu_check(coco_file_exists(file_name));
  mu_check(!coco_file_exists(nested));
  mu_check(!coco_directory_exists(file_name));
  mu_check(coco_create_unique_filename(nested, "data", ".txt") == "data-0001");

  /* Removing a file path is not removing a directory */
  mu_check(coco_remove_directory(file_name.c_str()) == 0);
  mu_check(coco_file_exists(file_name));

  /* Removes all contents */
  mu_check(coco_remove_directory(base.c_str()) == 0);
  mu_check(!coco_directory_exists(base));
  mu_check(coco_remove_directory(base.c_str()) == 0);

  /* Joining paths */
  mu_check(coco_join_path("", "x") == "x");
  mu_check(coco_join_path("a", "b", "c") == std::string("a") + (char)std::filesystem::path::preferred_separator + "b" +
                                               (char)std::filesystem::path::preferred_separator + "c");

}

/**
 * Tests the functions coco_double_max and coco_double_min.
 */
MU_TEST(test_coco_double_max_min) {

  double first_value = 5.0;
  double second_value = 6.0;

  double max_value, min_value;

  max_value = coco_double_max(first_value, second_value);
  mu_check(max_value == second_value);

  min_value = coco_double_min(first_value, second_value);
  mu_check(min_value == first_value);
}

/**
 * Tests the function coco_double_round.
 */
MU_TEST(test_coco_double_round) {

  double input_value = 5.4;

  double round_value = coco_double_round(input_value);
  mu_check(round_value == 5);

  input_value = 5.5;
  round_value = coco_double_round(input_value);
  mu_check(round_value == 6);

  input_value = 5;
  round_value = coco_double_round(input_value);
  mu_check(round_value == 5);

  input_value = 5 - 0.0000000000000001;
  round_value = coco_double_round(input_value);
  mu_check(round_value == 5);

  input_value = 5 + 0.0000000000000001;
  round_value = coco_double_round(input_value);
  mu_check(round_value == 5);
}

/**
 * Tests the function coco_double_round_up_with_precision.
 */
MU_TEST(test_coco_double_round_up_with_precision) {

  double input_value = 5.0024;
  double round_value = coco_double_round_up_with_precision(input_value, 0.00001);
  mu_check(coco_double_almost_equal(round_value, 5.00240, 1e-10));

  round_value = coco_double_round_up_with_precision(input_value, 0.0001);
  mu_check(coco_double_almost_equal(round_value, 5.0024, 1e-10));

  round_value = coco_double_round_up_with_precision(input_value, 0.001);
  mu_check(coco_double_almost_equal(round_value, 5.003, 1e-10));

  round_value = coco_double_round_up_with_precision(input_value, 0.01);
  mu_check(coco_double_almost_equal(round_value, 5.01, 1e-10));

  input_value = -5.0024;
  round_value = coco_double_round_up_with_precision(input_value, 0.00001);
  mu_check(coco_double_almost_equal(round_value, -5.00240, 1e-10));

  round_value = coco_double_round_up_with_precision(input_value, 0.0001);
  mu_check(coco_double_almost_equal(round_value, -5.0024, 1e-10));

  round_value = coco_double_round_up_with_precision(input_value, 0.001);
  mu_check(coco_double_almost_equal(round_value, -5.002, 1e-10));

  round_value = coco_double_round_up_with_precision(input_value, 0.01);
  mu_check(coco_double_almost_equal(round_value, -5, 1e-10));

  input_value = 0;
  round_value = coco_double_round_up_with_precision(input_value, 0.001);
  mu_check(coco_double_almost_equal(round_value, 0, 1e-10));

  input_value = -5.9999;
  round_value = coco_double_round_up_with_precision(input_value, 0.00001);
  mu_check(coco_double_almost_equal(round_value, -5.9999, 1e-10));

  round_value = coco_double_round_up_with_precision(input_value, 0.0001);
  mu_check(coco_double_almost_equal(round_value, -5.9999, 1e-10));

  round_value = coco_double_round_up_with_precision(input_value, 0.001);
  mu_check(coco_double_almost_equal(round_value, -5.999, 1e-10));

  input_value = 5.0001;
  round_value = coco_double_round_up_with_precision(input_value, 0.00001);
  mu_check(coco_double_almost_equal(round_value, 5.0001, 1e-10));

  round_value = coco_double_round_up_with_precision(input_value, 0.0001);
  mu_check(coco_double_almost_equal(round_value, 5.0001, 1e-10));

  round_value = coco_double_round_up_with_precision(input_value, 0.001);
  mu_check(coco_double_almost_equal(round_value, 5.001, 1e-10));

  input_value = 1e-8;
  round_value = coco_double_round_up_with_precision(input_value, 0.00001);
  mu_check(coco_double_almost_equal(round_value, 0.00001, 1e-10));

  round_value = coco_double_round_up_with_precision(input_value, 0.0001);
  mu_check(coco_double_almost_equal(round_value, 0.0001, 1e-10));

  round_value = coco_double_round_up_with_precision(input_value, 0.001);
  mu_check(coco_double_almost_equal(round_value, 0.001, 1e-10));

  input_value = 123;
  round_value = coco_double_round_up_with_precision(input_value, 20);
  mu_check(coco_double_almost_equal(round_value, 140, 1e-10));

  round_value = coco_double_round_up_with_precision(input_value, 120);
  mu_check(coco_double_almost_equal(round_value, 240, 1e-10));

  round_value = coco_double_round_up_with_precision(input_value, 123);
  mu_check(coco_double_almost_equal(round_value, 123, 1e-10));
}

/**
 * Tests the function coco_string_split.
 */
MU_TEST(test_coco_string_split) {

  mu_check(convert_to_string_with_newlines(coco_string_split("1-3,5-6,7,-3,15-", ',')) == "1-3\n5-6\n7\n-3\n15-\n");
  mu_check(convert_to_string_with_newlines(coco_string_split(",,,a,,b,c,d,e,,,,,f,", ',')) == "a\nb\nc\nd\ne\nf\n");
  mu_check(coco_string_split("", ',').empty());
  mu_check(coco_string_split(",,", ',').empty());
}

static std::string convert_to_string(const std::vector<size_t> &array) {

  std::string result;

  for (const size_t number : array) {
    result += coco_string_format("%lu,", (unsigned long)number);
  }
  return result + "0";
}

/**
 * Tests the function coco_string_parse_ranges.
 */
MU_TEST(test_coco_string_parse_ranges) {

  std::vector<size_t> result;

  mu_check(coco_string_parse_ranges("", 1, 2, "", 100).empty());
  mu_check(coco_string_parse_ranges("", 1, 2, "", 100).empty());
  mu_check(coco_string_parse_ranges("", 1, 2, "bla", 100).empty());
  mu_check(coco_string_parse_ranges("bla", 1, 2, "", 100).empty());
  mu_check(coco_string_parse_ranges("bla", 1, 2, "bla", 100).empty());

  result = coco_string_parse_ranges("1-3,5-6,7,-3,15-", 1, 20, "name", 100);
  mu_check(convert_to_string(result) == "1,2,3,5,6,7,1,2,3,15,16,17,18,19,20,0");

  result = coco_string_parse_ranges("1-3,5-6", 1, 2, "name", 100);
  mu_check(convert_to_string(result) == "1,2,0");

  result = coco_string_parse_ranges("1-3,5-6", 0, 20, "name", 100);
  mu_check(convert_to_string(result) == "1,2,3,5,6,0");

  result = coco_string_parse_ranges("1-", 1, 1200, "name", 1000);
  mu_check(result.size() == 1000);
  mu_check(result[999] == 1000);

  result = coco_string_parse_ranges("1-2,1-", 0, 0, "name", 1000);
  mu_check(convert_to_string(result) == "1,2,0");

  result = coco_string_parse_ranges("1-4,5-8", 2, 7, "name", 1000);
  mu_check(convert_to_string(result) == "2,3,4,5,6,7,0");

  mu_check(coco_string_parse_ranges("4-1", 1, 20, "name", 1000).empty());
  mu_check(coco_string_parse_ranges("-0", 1, 20, "name", 1000).empty());
  mu_check(coco_string_parse_ranges("-2,8-", 3, 7, "name", 1000).empty());

  result = coco_string_parse_ranges("1-8", 4, 8, "name", 1000);
  mu_check(convert_to_string(result) == "4,5,6,7,8,0");

  result = coco_string_parse_ranges(",7", 1, 20, "name", 1000);
  mu_check(convert_to_string(result) == "7,0");

  result = coco_string_parse_ranges("7,", 1, 20, "name", 1000);
  mu_check(convert_to_string(result) == "7,0");

  mu_check(coco_string_parse_ranges("-7-", 1, 20, "name", 1000).empty());
  mu_check(coco_string_parse_ranges("--7", 1, 20, "name", 1000).empty());
  mu_check(coco_string_parse_ranges("7--", 1, 20, "name", 1000).empty());

  result = coco_string_parse_ranges("7-,-", 5, 8, "name", 1000);
  mu_check(convert_to_string(result) == "7,8,5,6,7,8,0");

  result = coco_string_parse_ranges("7-,,5", 1, 8, "name", 1000);
  mu_check(convert_to_string(result) == "7,8,5,0");

  mu_check(coco_string_parse_ranges(",,", 1, 20, "name", 1000).empty());
  mu_check(coco_string_parse_ranges("1-8", 5, 2, "name", 1000).empty());
}

/**
 * Tests the function coco_option_keys.
 */
MU_TEST(test_coco_option_keys) {

  coco_option_keys_t option_keys;

  mu_check(coco_option_keys("").empty());
  mu_check(coco_option_keys("").empty());

  option_keys = coco_option_keys("key");
  mu_check(option_keys.size() == 1);
  mu_check(option_keys[0] == "key");

  option_keys = coco_option_keys("key: ");
  mu_check(option_keys.size() == 1);

  option_keys = coco_option_keys("key1 key2: ");
  mu_check(option_keys.size() == 1);

  option_keys = coco_option_keys("key1: value1 key2");
  /* In this case we would rather have detected two keys, but this should also trigger a warning,
   * which is OK. */
  mu_check(option_keys.size() == 1);

  option_keys = coco_option_keys("key1: value1 key2: value2");
  mu_check(option_keys.size() == 2);
  mu_check(option_keys[0] == "key1");
  mu_check(option_keys[1] == "key2");

  option_keys = coco_option_keys("key: \"A multi-word value\"");
  mu_check(option_keys.size() == 1);

  option_keys = coco_option_keys("key: \"A multi-word value with a colon (:)\"");
  mu_check(option_keys.size() == 1);

  option_keys = coco_option_keys("\"key: value\"");
  mu_check(option_keys.size() == 1);

  option_keys = coco_option_keys_get_redundant({"key1", "key2"}, coco_option_keys("key1: 1 key3: 3 key2: 2"));
  mu_check(option_keys.size() == 1);
  mu_check(option_keys[0] == "key3");
  mu_check(coco_option_keys_get_output_string(option_keys, "Info:\n") == "Info:\n key3\n");
  mu_check(coco_option_keys_get_output_string(coco_option_keys_t(), "Info:\n").empty());
}

/**
 * Tests the functions coco_options_read_string and coco_options_read_values.
 */
MU_TEST(test_coco_options_read) {

  std::string value = "default";

  mu_check(coco_options_read_string("name: value other: 2", "name", value) == 1);
  mu_check(value == "value");
  mu_check(coco_options_read_string("name: \"multi word value\" other: 2", "name", value) == 1);
  mu_check(value == "multi word value");
  value = "default";
  mu_check(coco_options_read_string("other: 2", "name", value) == 0);
  mu_check(coco_options_read_string("name:", "name", value) == 0);
  mu_check(coco_options_read_string("", "name", value) == 0);
  mu_check(value == "default");

  mu_check(coco_options_read_values("indices: 1, 2-4 other: 3", "indices", value) == 5);
  mu_check(value == "1,2-4");
  mu_check(coco_options_read_values("other: 3", "indices", value) == 0);
  mu_check(value.empty());
}

/**
 * Tests the function coco_is_nan.
 */
MU_TEST(test_coco_is_nan) {

  mu_check(coco_is_nan(NAN));
  mu_check(!coco_is_nan(1e100));
}

/**
 * Tests the function coco_is_inf.
 */
MU_TEST(test_coco_is_inf) {

  mu_check(coco_is_inf(INFINITY));
  mu_check(coco_is_inf(-INFINITY));
  mu_check(coco_is_inf(2*INFINITY));
  mu_check(coco_is_inf(-2*INFINITY));
  mu_check(!coco_is_inf(NAN));
}

/**
 * Tests coco_matrix_t::solve().
 */
MU_TEST(test_coco_matrix_solve) {
  /* Needs a row exchange because of the zero in the upper left corner */
  coco_matrix_t M(3, 3);
  const double values[3][3] = {{0., 2., 1.}, {1., 1., 0.}, {4., 0., 2.}};
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      M(i, j) = values[i][j];
  const double x_expected[3] = {1., -2., 3.};
  double y[3], x[3];
  M.multiply(x_expected, y);
  mu_check(M.solve(y, x));
  for (size_t i = 0; i < 3; ++i)
    mu_check(fabs(x[i] - x_expected[i]) < 1e-14);

  /* A singular matrix */
  for (size_t j = 0; j < 3; ++j)
    M(2, j) = M(0, j) + M(1, j);
  mu_check(!M.solve(y, x));

  /* A non-square matrix */
  mu_check(!coco_matrix_t(2, 3).solve(y, x));
}

/**
 * Tests coco_matrix_t::is_orthogonal().
 */
MU_TEST(test_coco_matrix_is_orthogonal) {

  /* The identity matrix */
  coco_matrix_t M(2, 2);
  M(0, 0) = 1.;
  M(1, 1) = 1.;
  mu_check(M.is_orthogonal());

  /* A square matrix that is not orthogonal */
  M(0, 1) = 0.01;
  mu_check(!M.is_orthogonal());

  /* Close enough to the identity matrix (tests the use of coco_double_almost_equal()) */
  M(0, 0) = 1. + 1e-10;
  M(0, 1) = 0.;
  M(1, 1) = 1. + 1e-10;
  mu_check(M.is_orthogonal());

  /* A rotation by 90 degrees is orthogonal, a scaled one is not */
  M(0, 0) = 0.;
  M(0, 1) = -1.;
  M(1, 0) = 1.;
  M(1, 1) = 0.;
  mu_check(M.is_orthogonal());
  M(1, 0) = 2.;
  mu_check(!M.is_orthogonal());

  /* A non-square matrix is never orthogonal, even with orthonormal rows */
  coco_matrix_t rows(2, 3);
  rows(0, 0) = 1.;
  rows(1, 1) = 1.;
  mu_check(!rows.is_orthogonal());

  /* The rotations used by the bbob functions are orthogonal */
  mu_check(bbob2009_compute_rotation(42, 10).is_orthogonal());
}

/**
 * Tests coco_matrix_t.
 */
MU_TEST(test_coco_matrix) {
  coco_matrix_t matrix(2, 3);
  const double x[3] = {1.0, 2.0, 3.0};
  double y[2];

  mu_check(matrix.number_of_rows() == 2 && matrix.number_of_columns() == 3);
  mu_check(matrix(1, 2) == 0.0);
  matrix(0, 0) = 1.0;
  matrix(1, 2) = 5.0;
  mu_check(matrix.row(1)[2] == 5.0);
  mu_check(matrix.row(1).data() == matrix.row(0).data() + 3); /* row-major */

  matrix.multiply(x, y);
  mu_check(y[0] == 1.0 && y[1] == 15.0);

  /* Copies are independent */
  coco_matrix_t copy = matrix;
  copy(1, 2) = 7.0;
  mu_check(matrix(1, 2) == 5.0);
}

/**
 * Tests coco_block_diagonal_matrix_t.
 */
MU_TEST(test_coco_block_diagonal_matrix) {
  const size_t block_sizes[2] = {1, 2};
  coco_block_diagonal_matrix_t matrix(block_sizes);
  const double x[3] = {1.0, 2.0, 3.0};
  double y[3];

  mu_check(matrix.dimension() == 3);
  mu_check(matrix.block_size(0) == 1 && matrix.block_size(2) == 2);
  mu_check(matrix.first_column(0) == 0 && matrix.first_column(1) == 1 && matrix.first_column(2) == 1);
  matrix.in_block(0, 0) = 2.0;
  matrix.in_block(1, 1) = 3.0;
  matrix.in_block(2, 0) = 4.0;
  mu_check(matrix(1, 2) == 3.0 && matrix(2, 1) == 4.0 && matrix(0, 1) == 0.0 && matrix(1, 0) == 0.0);

  matrix.multiply(x, y);
  mu_check(y[0] == 2.0 && y[1] == 9.0 && y[2] == 8.0);
}

/**
 * Tests coco_string_format.
 */
MU_TEST(test_coco_string_format) {
  std::string long_string(1000, 'x');

  mu_check(coco_string_format("%s_d%02lu", "sphere", 5ul) == "sphere_d05");
  mu_check(coco_string_format("%s", "") == "");
  /* No length limit */
  mu_check(coco_string_format("%s!", long_string.c_str()).size() == 1001);
}

/**
 * Run all tests in this file.
 */
int main(void) {
  MU_RUN_TEST(test_coco_matrix);
  MU_RUN_TEST(test_coco_block_diagonal_matrix);
  MU_RUN_TEST(test_coco_string_format);
  MU_RUN_TEST(test_coco_is_nan);
  MU_RUN_TEST(test_coco_is_inf);
  MU_RUN_TEST(test_coco_set_log_level);
  MU_RUN_TEST(test_coco_double_max_min);
  MU_RUN_TEST(test_coco_double_round);
  MU_RUN_TEST(test_coco_double_round_up_with_precision);
  MU_RUN_TEST(test_coco_string_split);
  MU_RUN_TEST(test_coco_option_keys);
  MU_RUN_TEST(test_coco_string_parse_ranges);
  MU_RUN_TEST(test_coco_options_read);
  MU_RUN_TEST(test_coco_create_remove_directory);
  MU_RUN_TEST(test_coco_directories_and_files);
  MU_RUN_TEST(test_coco_matrix_is_orthogonal);
  MU_RUN_TEST(test_coco_matrix_solve);

  MU_REPORT();

  return MU_EXIT_CODE;
}
