#ifndef SCPH_LANDSCAPE_KERNEL_H
#define SCPH_LANDSCAPE_KERNEL_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum scph_landscape_status_v1 {
  SCPH_LANDSCAPE_OK_V1 = 0,
  SCPH_LANDSCAPE_NULL_OUTPUT_V1 = 1,
  SCPH_LANDSCAPE_NULL_INPUT_V1 = 2,
  SCPH_LANDSCAPE_INVALID_INPUT_V1 = 3,
  SCPH_LANDSCAPE_PANIC_V1 = 4,
  SCPH_LANDSCAPE_NONFINITE_OUTPUT_V1 = 5,
  SCPH_LANDSCAPE_UNSUPPORTED_DIMENSION_V1 = 6
};

struct scph_landscape_result_v1 {
  double squared_distance;
  uint64_t active_levels;
  uint64_t event_segments;
  uint64_t first_finite_intervals;
  uint64_t second_finite_intervals;
  uint32_t engine_version;
  int32_t status;
};

int32_t scph_landscape_l2_v1(
    const double *first_births,
    const double *first_deaths,
    size_t first_len,
    const double *second_births,
    const double *second_deaths,
    size_t second_len,
    uint32_t dimension,
    struct scph_landscape_result_v1 *output);

#ifdef __cplusplus
}
#endif

#endif
