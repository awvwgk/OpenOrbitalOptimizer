/*
 Copyright (C) 2023- Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

// Exercises the Aufbau cleanup on a Fermi level that is degenerate
// across more than two orbitals, which is where the cleanup has to
// model the occupation-orbital coupling over a polytope rather than a
// line.
//
// The atomic driver cannot reach this. Its blocks hold whole l shells
// with capacities 2(2l+1), so a degenerate 4s/3d pair admits only the
// two skeleton fillings and the polytope is a line. A block of
// singly-occupiable orbitals is what produces more: three of them
// sharing two particles admit (1,1,0), (1,0,1) and (0,1,1), so the
// polytope is a triangle. Nothing about that is exotic -- it is what a
// caller gets whenever a block resolves individual orbitals rather
// than shells.

#include <openorbitaloptimizer/scfsolver.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using OpenOrbitalOptimizer::SCFSolver;
using OpenOrbitalOptimizer::Matrix;
using OpenOrbitalOptimizer::Vector;
using OpenOrbitalOptimizer::IndexVector;
using OpenOrbitalOptimizer::FockBuilder;
using OpenOrbitalOptimizer::FockBuilderReturn;
using OpenOrbitalOptimizer::DensityMatrix;
using OpenOrbitalOptimizer::FockMatrix;

namespace {

int failures = 0;

#define REQUIRE(cond)                                                   \
  do {                                                                  \
    if (!(cond)) {                                                      \
      std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);       \
      ++failures;                                                       \
    }                                                                   \
  } while (0)

// A model whose Fermi level is degenerate across three orbitals, with
// genuine orbital freedom underneath it.
//
//   E = tr(h D) + (u/2) sum_p D_pp^2
//
// h is the hopping matrix of the complete graph on n_orbital sites,
// which puts one level well below a cluster of n_orbital-1 that are
// degenerate by permutation symmetry, plus a small diagonal splitting
// so that cluster is only degenerate to within the solver's window
// rather than exactly. That matters: an exactly symmetric model has
// its orbitals fixed by symmetry, F commuting with h, so there would
// be nothing for the relaxation to do and the occupations and orbitals
// would not be coupled at all. Splitting the level slightly leaves the
// cluster intact -- the window is 1e-2 and the splitting 1e-3 -- while
// giving the orbitals something to relax.
//
// The on-site repulsion is what makes the energy non-linear in the
// occupations and, being written in the site basis, does not commute
// with h, so the orbitals genuinely respond to an occupation change.
SCFSolver<double, double> make_degenerate_solver(int n_orbital,
                                                 double n_particles,
                                                 double u) {
  IndexVector blocks_per_particle(1);
  blocks_per_particle(0) = 1;
  Vector<double> maxocc(1);
  maxocc(0) = 1.0;                    // one particle per orbital
  Vector<double> nparticles(1);
  nparticles(0) = n_particles;

  Matrix<double> h = Matrix<double>::Zero(n_orbital, n_orbital);
  for (int p = 0; p < n_orbital; ++p)
    for (int q = 0; q < n_orbital; ++q)
      if (p != q) h(p, q) = -1.0;
  // Split the degenerate cluster well inside the degeneracy window.
  for (int p = 1; p < n_orbital; ++p)
    h(p, p) = 1e-3 * p;

  FockBuilder<double, double> fb =
    [n_orbital, h, u](const DensityMatrix<double, double> & dm) {
      const auto & C = dm.first[0];
      const auto & n = dm.second[0];

      Matrix<double> D = Matrix<double>::Zero(n_orbital, n_orbital);
      for (int k = 0; k < n_orbital; ++k)
        D += n(k) * C.col(k) * C.col(k).transpose();

      double energy = (h * D).trace();
      Matrix<double> F = h;
      for (int p = 0; p < n_orbital; ++p) {
        energy += 0.5 * u * D(p, p) * D(p, p);
        F(p, p) += u * D(p, p);
      }
      return FockBuilderReturn<double, double>{energy, {F}};
    };

  return SCFSolver<double, double>(blocks_per_particle, maxocc, nparticles,
                                   fb, {"level"});
}

// Largest occupation above the Fermi level or missing from below it,
// read off the solver's own diagnostic.
double aufbau_error(const SCFSolver<double, double> & s) {
  return s.get_real("aufbau_error");
}

}  // namespace

int main() {
  // Three degenerate orbitals sharing two particles: the skeleton set
  // is {(1,1,0), (1,0,1), (0,1,1)} and the polytope a triangle, so the
  // cleanup takes the multi-dimensional branch.
  {
    auto s = make_degenerate_solver(/*n_orbital=*/4, /*n_particles=*/2.0,
                                    /*u=*/0.3);
    FockMatrix<double> guess(1);
    guess[0] = Matrix<double>::Zero(4, 4);
    for (int p = 0; p < 4; ++p)
      for (int q = 0; q < 4; ++q)
        if (p != q) guess[0](p, q) = -1.0;
    s.initialize_with_fock(guess);
    s.set(std::string("verbosity"), 0);
    s.set(std::string("methods"), std::string("ODA + LBFGS"));
    s.set(std::string("maximum_iterations"), 200);
    s.run();

    REQUIRE(s.converged());
    // Particle number is an invariant of every step; nothing may
    // discard any of it.
    REQUIRE(s.particle_number_error() <= 1e-10);
    // And the reported occupations have to be Aufbau, which is the
    // whole point of the cleanup.
    REQUIRE(aufbau_error(s) <= s.get_real("aufbau_convergence_threshold"));

    const auto occupations = s.get_orbital_occupations();
    std::printf("  converged occupations:");
    for (int k = 0; k < occupations[0].size(); ++k)
      std::printf(" %.6f", occupations[0](k));
    std::printf("   E = %.10f  aufbau_error = %.3e\n",
                s.get_energy(), aufbau_error(s));
    double total = 0.0;
    for (int k = 0; k < occupations[0].size(); ++k) {
      const double n = occupations[0](k);
      total += n;
      // Occupations must stay inside [0, max_occ]; the polytope search
      // has no business leaving the simplex.
      REQUIRE(n >= -1e-10);
      REQUIRE(n <= 1.0 + 1e-10);
    }
    REQUIRE(std::abs(total - 2.0) <= 1e-10);
  }

  // Four degenerate orbitals sharing two particles, a larger polytope
  // over the same machinery.
  {
    auto s = make_degenerate_solver(/*n_orbital=*/5, /*n_particles=*/2.0,
                                    /*u=*/0.3);
    FockMatrix<double> guess(1);
    guess[0] = Matrix<double>::Zero(5, 5);
    for (int p = 0; p < 5; ++p)
      for (int q = 0; q < 5; ++q)
        if (p != q) guess[0](p, q) = -1.0;
    s.initialize_with_fock(guess);
    s.set(std::string("verbosity"), 0);
    s.set(std::string("methods"), std::string("ODA + LBFGS"));
    s.set(std::string("maximum_iterations"), 200);
    s.run();

    REQUIRE(s.converged());
    REQUIRE(s.particle_number_error() <= 1e-10);
    REQUIRE(aufbau_error(s) <= s.get_real("aufbau_convergence_threshold"));

    const auto occupations = s.get_orbital_occupations();
    std::printf("  converged occupations:");
    for (int k = 0; k < occupations[0].size(); ++k)
      std::printf(" %.6f", occupations[0](k));
    std::printf("   E = %.10f  aufbau_error = %.3e\n",
                s.get_energy(), aufbau_error(s));
    double total = 0.0;
    for (int k = 0; k < occupations[0].size(); ++k) {
      const double n = occupations[0](k);
      total += n;
      REQUIRE(n >= -1e-10);
      REQUIRE(n <= 1.0 + 1e-10);
    }
    REQUIRE(std::abs(total - 2.0) <= 1e-10);
  }

  if (failures) {
    std::printf("%d failure(s)\n", failures);
    return 1;
  }
  std::printf("All degenerate-occupation checks passed\n");
  return 0;
}
