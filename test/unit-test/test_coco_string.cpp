#include "minunit.h"

#include "coco.cpp"

/**
 * Tests coco_string-related functions.
 */
MU_TEST(test_coco_string_trim) {

	const char *sample_strings[] = {
			"nothing to trim",
			"    trim the front",
			"trim the back     ",
			" trim one char front and back ",
			" trim one char front",
			"trim one char back ",
			"                   ",
			" ",
			"a",
			"",
			NULL };
	const char *result_strings[] = {
			"nothing to trim",
			"trim the front",
			"trim the back",
			"trim one char front and back",
			"trim one char front",
			"trim one char back",
			"",
			"",
			"a",
			"",
			NULL };
	int index;

	for (index = 0; sample_strings[index] != NULL; ++index) {
		mu_check(coco_string_trim(sample_strings[index]) == result_strings[index]);
	}
}

/**
 * Run all tests in this file.
 */
int main(void) {
  MU_RUN_TEST(test_coco_string_trim);

  MU_REPORT();

  return MU_EXIT_CODE;
}
