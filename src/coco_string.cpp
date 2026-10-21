/**
 * @file coco_string.cpp
 * @brief Definitions of functions that manipulate strings.
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>

#include <algorithm>
#include <string>
#include <string_view>
#include <vector>

#include "coco.h"

/**
 * @brief Wraps a string that comes in through the C interface; NULL becomes an empty string.
 */
static std::string_view coco_c_string_view(const char *string) {
  return (string == nullptr) ? std::string_view() : std::string_view(string);
}

/**
 * @brief Formats a string (as in vsprintf) and returns it as std::string.
 */
static std::string coco_string_vformat(const char *format, va_list args) {
  va_list args_copy;
  int length;
  std::string result;

  va_copy(args_copy, args);
  length = vsnprintf(nullptr, 0, format, args_copy);
  va_end(args_copy);
  if (length < 0)
    coco_error("coco_string_vformat(): vsnprintf failed on '%s'", format);
  result.resize((size_t)length + 1);
  vsnprintf(&result[0], result.size(), format, args);
  result.resize((size_t)length);
  return result;
}

/**
 * @brief Formats a string (as in sprintf) and returns it as std::string.
 */
COCO_FORMAT(1, 2) static std::string coco_string_format(const char *format, ...) {
  va_list args;
  std::string result;

  va_start(args, format);
  result = coco_string_vformat(format, args);
  va_end(args);
  return result;
}

/**
 * Optional arguments are used like in sprintf. The caller frees the result with coco_free_memory().
 */
char *coco_strdupf(const char *str, ...) {
  va_list args;

  va_start(args, str);
  const std::string string = coco_string_vformat(str, args);
  va_end(args);

  char *duplicate = (char *)coco_allocate_memory(string.size() + 1);
  memcpy(duplicate, string.c_str(), string.size() + 1);
  return duplicate;
}

/**
 * @brief Returns the first index where seq occurs in base and -1 if it doesn't.
 */
static long coco_strfind(std::string_view base, std::string_view seq) {
  const size_t index = base.find(seq);
  return (index == std::string_view::npos) ? -1 : (long)index;
}

/**
 * @brief Splits a string based on the given delimiter.
 *
 * Empty substrings (caused by consecutive delimiters or delimiters at the start or end) are skipped.
 */
static std::vector<std::string> coco_string_split(std::string_view string, const char delimiter) {
  std::vector<std::string> result;
  size_t start = 0;

  while (start < string.size()) {
    size_t stop = string.find(delimiter, start);
    if (stop == std::string_view::npos)
      stop = string.size();
    if (stop > start)
      result.emplace_back(string.substr(start, stop - start));
    start = stop + 1;
  }
  return result;
}

/**
 * @brief Returns a copy of the string with the characters between from and to removed.
 *
 * If you wish to remove characters from the beginning of the string, set from to "".
 * If you wish to remove characters until the end of the string, set to to "".
 */
static std::string coco_remove_from_string(const std::string &string, std::string_view from, std::string_view to) {

  const size_t start = from.empty() ? 0 : string.find(from);
  const size_t stop = to.empty() ? string.size() : string.find(to);

  if ((start == std::string::npos) || (stop == std::string::npos) || (stop < start)) {
    coco_error("coco_remove_from_string(): failed to remove characters between %.*s and %.*s from string %s",
               (int)from.size(), from.data(), (int)to.size(), to.data(), string.c_str());
    return std::string(); /* Never reached */
  }
  return string.substr(0, start) + string.substr(stop);
}

/**
 * @brief Returns the numbers defined by the ranges.
 *
 * Reads ranges from a string of positive ranges separated by commas. For example: "-3,5-6,8-". Returns the
 * numbers that are defined by the ranges if min and max are used as their extremes. If the ranges with open
 * beginning/end are not allowed, use 0 as min/max. A maximum of max_count values is returned. If there is a
 * problem with one of the ranges, the parsing stops and the current result is returned. An empty result
 * means that no numbers could be parsed.
 */
static std::vector<size_t> coco_string_parse_ranges(std::string_view string, const size_t min, const size_t max,
                                                    const char *name, const size_t max_count) {

  std::vector<size_t> result;
  size_t j, num[2];

  /* Check for empty string */
  if (string.empty()) {
    coco_warning("coco_string_parse_ranges(): cannot parse empty ranges");
    return result;
  }

  /* Check for disallowed characters */
  for (const char c : string) {
    if ((c != '-') && (c != ',') && !isdigit((unsigned char)c)) {
      coco_warning("coco_string_parse_ranges(): problem parsing '%.*s' - cannot parse ranges with '%c'",
                   (int)string.size(), string.data(), c);
      return result;
    }
  }

  /* Check for incorrect boundaries */
  if ((max > 0) && (min > max)) {
    coco_warning("coco_string_parse_ranges(): incorrect boundaries");
    return result;
  }

  /* Split string to ranges w.r.t commas and go over them */
  for (const std::string &range : coco_string_split(string, ',')) {

    const size_t count = (size_t)std::count(range.begin(), range.end(), '-');

    if (count > 1) {
      /* Incorrect number of '-' */
      coco_warning("coco_string_parse_ranges(): problem parsing '%.*s' - too many '-'s", (int)string.size(),
                   string.data());
      return result;
    } else if (count == 0) {
      /* Range is in the format: n (no range) */
      num[0] = (size_t)strtol(range.c_str(), nullptr, 10);
      num[1] = num[0];
    } else {
      /* Range is in one of the following formats: n-m / -n / n- / - */
      const std::vector<std::string> numbers = coco_string_split(range, '-');
      for (j = 0; j < numbers.size(); j++) {
        assert(j < 2);
        num[j] = (size_t)strtol(numbers[j].c_str(), nullptr, 10);
      }

      if (numbers.empty()) {
        /* Range is in the format - (open ends) */
        if ((min == 0) || (max == 0)) {
          coco_warning("coco_string_parse_ranges(): '%s' ranges cannot have an open ends; some ranges ignored", name);
          return result;
        }
        num[0] = min;
        num[1] = max;
      } else if (numbers.size() == 1) {
        if (range[0] == '-') {
          /* Range is in the format -n */
          if (min == 0) {
            coco_warning("coco_string_parse_ranges(): '%s' ranges cannot have an open beginning; some ranges ignored",
                         name);
            return result;
          }
          num[1] = num[0];
          num[0] = min;
        } else {
          /* Range is in the format n- */
          if (max == 0) {
            coco_warning("coco_string_parse_ranges(): '%s' ranges cannot have an open end; some ranges ignored", name);
            return result;
          }
          num[1] = max;
        }
      }
      /* if (numbers.size() == 2), range is in the format n-m and there is nothing to do */
    }

    /* Make sure the boundaries are taken into account */
    if ((min > 0) && (num[0] < min)) {
      num[0] = min;
      coco_warning("coco_string_parse_ranges(): '%s' ranges adjusted to be >= %lu", name, (unsigned long)min);
    }
    if ((max > 0) && (num[1] > max)) {
      num[1] = max;
      coco_warning("coco_string_parse_ranges(): '%s' ranges adjusted to be <= %lu", name, (unsigned long)max);
    }
    if (num[0] > num[1]) {
      coco_warning("coco_string_parse_ranges(): '%s' ranges not within boundaries; some ranges ignored", name);
      return result;
    }

    /* Write in result */
    for (j = num[0]; (j <= num[1]) && (result.size() < max_count); j++) {
      result.push_back(j);
    }
  }

  return result;
}

/**
 * @brief Returns the given string without leading and trailing white space.
 */
static std::string coco_string_trim(std::string_view string) {
  size_t start = 0, stop = string.size();

  while ((start < stop) && isspace((unsigned char)string[start]))
    start++;
  while ((stop > start) && isspace((unsigned char)string[stop - 1]))
    stop--;
  return std::string(string.substr(start, stop - start));
}
