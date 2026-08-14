#include "scph_landscape_kernel.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>

int main(void) {
  const double first_birth[] = {0.0};
  const double first_death[] = {2.0};
  const double second_birth[] = {0.25};
  const double second_death[] = {2.25};
  struct scph_landscape_result_v1 result = {0};

  int32_t status = scph_landscape_l2_v1(
      first_birth, first_death, 1, second_birth, second_death, 1, 0,
      &result);
  assert(status == SCPH_LANDSCAPE_OK_V1);
  assert(result.status == SCPH_LANDSCAPE_OK_V1);
  assert(result.engine_version == 1);
  assert(fabs(result.squared_distance - 7.0 / 64.0) <= 1e-15);
  assert(result.active_levels == 1);

  status = scph_landscape_l2_v1(
      first_birth, first_birth, 1, NULL, NULL, 0, 0, &result);
  assert(status == SCPH_LANDSCAPE_INVALID_INPUT_V1);
  assert(result.status == SCPH_LANDSCAPE_INVALID_INPUT_V1);

  status = scph_landscape_l2_v1(
      first_birth, first_death, 1, NULL, NULL, 0, 2, &result);
  assert(status == SCPH_LANDSCAPE_UNSUPPORTED_DIMENSION_V1);
  assert(result.status == SCPH_LANDSCAPE_UNSUPPORTED_DIMENSION_V1);

  status = scph_landscape_l2_v1(
      first_birth, first_death, 1, NULL, NULL, 0, 0, NULL);
  assert(status == SCPH_LANDSCAPE_NULL_OUTPUT_V1);

  status = scph_landscape_l2_v1(
      NULL, NULL, 1, NULL, NULL, 0, 0, &result);
  assert(status == SCPH_LANDSCAPE_NULL_INPUT_V1);
  assert(result.status == SCPH_LANDSCAPE_NULL_INPUT_V1);
  return 0;
}
