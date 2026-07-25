#include <openorbitaloptimizer/quad_support.hpp>
#include <openorbitaloptimizer/scfsolver.hpp>
#include "settings_api_check.hpp"

template class OpenOrbitalOptimizer::SCFSolver<std::complex<_Float128>, _Float128>;

int main(void) {
  // Force instantiation of the settings-facade member templates,
  // which explicit class instantiation does not cover.
  auto * p = &OpenOrbitalOptimizer::SettingsApiCheck::check<std::complex<_Float128>,_Float128>;
  (void) p;
  return 0;
}
