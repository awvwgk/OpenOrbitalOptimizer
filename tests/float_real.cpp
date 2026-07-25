#include <openorbitaloptimizer/scfsolver.hpp>
#include "settings_api_check.hpp"

template class OpenOrbitalOptimizer::SCFSolver<float, float>;

int main(void) {
  // Force instantiation of the settings-facade member templates,
  // which explicit class instantiation does not cover.
  auto * p = &OpenOrbitalOptimizer::SettingsApiCheck::check<float,float>;
  (void) p;
  return 0;
}
