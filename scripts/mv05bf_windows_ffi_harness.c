#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include "scph_landscape_kernel.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

typedef int32_t(__cdecl *scph_landscape_l2_fn)(
    const double *, const double *, size_t, const double *, const double *,
    size_t, uint32_t, struct scph_landscape_result_v1 *);

int main(int argc, char **argv) {
  HMODULE library;
  scph_landscape_l2_fn landscape;
  FARPROC address;
  const double first_birth[] = {0.0};
  const double first_death[] = {2.0};
  const double second_birth[] = {0.25};
  const double second_death[] = {2.25};
  struct scph_landscape_result_v1 result = {0};
  int32_t status;

  assert(argc == 2);
  library = LoadLibraryA(argv[1]);
  assert(library != NULL);
  address = GetProcAddress(library, "scph_landscape_l2_v1");
  assert(address != NULL);
  assert(sizeof(address) == sizeof(landscape));
  memcpy(&landscape, &address, sizeof(landscape));

  status = landscape(first_birth, first_death, 1, second_birth, second_death,
                     1, 0, &result);
  assert(status == SCPH_LANDSCAPE_OK_V1);
  assert(result.status == SCPH_LANDSCAPE_OK_V1);
  assert(result.engine_version == 1);
  assert(fabs(result.squared_distance - 7.0 / 64.0) <= 1e-15);

  status = landscape(first_birth, first_birth, 1, NULL, NULL, 0, 0, &result);
  assert(status == SCPH_LANDSCAPE_INVALID_INPUT_V1);
  status = landscape(first_birth, first_death, 1, NULL, NULL, 0, 2, &result);
  assert(status == SCPH_LANDSCAPE_UNSUPPORTED_DIMENSION_V1);
  status = landscape(first_birth, first_death, 1, NULL, NULL, 0, 0, NULL);
  assert(status == SCPH_LANDSCAPE_NULL_OUTPUT_V1);

  assert(FreeLibrary(library));
  return 0;
}
