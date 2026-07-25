#include <openorbitaloptimizer/scfsolver.hpp>
#include "settings_api_check.hpp"

// Instantiate the SCF solver class
template class OpenOrbitalOptimizer::SCFSolver<std::complex<double>, double>;

int main(void) {
  // Force instantiation of the settings-facade member templates,
  // which explicit class instantiation does not cover.
  auto * p = &OpenOrbitalOptimizer::SettingsApiCheck::check<std::complex<double>,double>;
  (void) p;
  return 0;
}
