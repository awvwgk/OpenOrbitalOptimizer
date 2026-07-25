/*
 Copyright (C) 2023- Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

// Exercises the SCFSolver settings façade: options() catalog,
// set(key, value) / get_*(key) round-trip on every entry, and
// invalid-key rejection. Runs at ctest time so a stale catalog or
// dispatch mismatch trips CI before it reaches downstream callers.

#include <openorbitaloptimizer/scfsolver.hpp>

#include <cstdio>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <stdexcept>
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

SCFSolver<double, double> make_solver() {
  IndexVector blocks_per_particle(1);
  blocks_per_particle(0) = 1;
  Vector<double> maxocc(1);
  maxocc(0) = 2.0;
  Vector<double> nparticles(1);
  nparticles(0) = 2.0;
  FockBuilder<double, double> fb =
    [](const DensityMatrix<double, double> &) {
      FockMatrix<double> F(1);
      F[0] = Matrix<double>::Identity(2, 2) * -1.0;
      return FockBuilderReturn<double, double>{0.0, F};
    };
  return SCFSolver<double, double>(blocks_per_particle, maxocc, nparticles,
                                   fb, {"s"});
}

}  // namespace

int main() {
  auto solver = make_solver();

  // Seed the history so read-only diagnostics that walk it
  // (converged, DIIS-error-based ones) have data to read.
  FockMatrix<double> guess(1);
  guess[0] = Matrix<double>::Identity(2, 2) * -1.0;
  solver.initialize_with_fock(guess);

  const auto & catalog = SCFSolver<double, double>::options();
  REQUIRE(!catalog.empty());

  // Round-trip every writable entry and read every read-only entry.
  for (const auto & opt : catalog) {
    const std::string key = opt.key;
    const std::string type = opt.type;
    if (opt.writable) {
      if (type == "real") {
        double v = 0.1234567;
        solver.set(key, v);
        REQUIRE(solver.get_real(key) == v);
      } else if (type == "int") {
        // frozen_occupations rides on int but coerces to bool; 1
        // round-trips through both cleanly for every int entry.
        int v = 1;
        solver.set(key, v);
        REQUIRE(solver.get_int(key) == v);
      } else if (type == "string") {
        // Every catalog string is validated on set(); pick an input
        // we know each accepts. `methods` is stored in canonical
        // uppercase, so the round-trip echoes back "ODA".
        std::string v_in, v_expect;
        if (key == "error_norm")   { v_in = "rms";  v_expect = "rms"; }
        else if (key == "methods") { v_in = "oda";  v_expect = "ODA"; }
        else                       { v_in = "";     v_expect = ""; }
        solver.set(key, v_in);
        REQUIRE(solver.get_string(key) == v_expect);
      } else {
        std::printf("FAIL: unknown type '%s' in catalog for '%s'\n",
                    type.c_str(), key.c_str());
        ++failures;
      }
    } else {
      // Read-only: just make sure the getter doesn't throw.
      try {
        if      (type == "real")   (void) solver.get_real(key);
        else if (type == "int")    (void) solver.get_int(key);
        else if (type == "string") (void) solver.get_string(key);
      } catch (const std::exception & e) {
        std::printf("FAIL: read-only get for '%s' threw: %s\n",
                    key.c_str(), e.what());
        ++failures;
      }
    }
  }

  // Unknown key rejection.
  bool threw = false;
  try { solver.set(std::string("no_such_setting"), 1.0); }
  catch (const std::invalid_argument &) { threw = true; }
  REQUIRE(threw);

  threw = false;
  try { (void) solver.get_int("no_such_setting"); }
  catch (const std::invalid_argument &) { threw = true; }
  REQUIRE(threw);

  // Wrong-type rejection: convergence_threshold is real, not int.
  threw = false;
  try { (void) solver.get_int("convergence_threshold"); }
  catch (const std::invalid_argument &) { threw = true; }
  REQUIRE(threw);

  // print_settings covers every catalog entry; every key appears in the
  // output, and print_citation names the paper.
  {
    std::ostringstream oss;
    solver.print_settings(oss);
    std::string dump = oss.str();
    for (const auto & o : catalog) {
      if (dump.find(o.key) == std::string::npos) {
        std::printf("FAIL: print_settings missed key '%s'\n", o.key);
        ++failures;
      }
    }
  }
  {
    std::ostringstream oss;
    SCFSolver<double, double>::print_citation(oss);
    std::string cite = oss.str();
    REQUIRE(cite.find("Lehtola") != std::string::npos);
    REQUIRE(cite.find("10.1021/acs.jpca.5c02110") != std::string::npos);
    REQUIRE(cite.find("J. Phys. Chem. A") != std::string::npos);
  }

  // Logger callback catches messages instead of stdout, and honours
  // the verbosity gate.
  {
    auto logger_solver = make_solver();
    FockMatrix<double> guess(1);
    guess[0] = Matrix<double>::Identity(2, 2) * -1.0;
    logger_solver.initialize_with_fock(guess);

    struct Record {
      int level;
      std::string msg;
    };
    std::vector<Record> captured;
    logger_solver.logger(
        [&](int level, const std::string & msg) {
          captured.push_back({level, msg});
        });
    REQUIRE(logger_solver.has_logger());

    // A run() at the default verbosity emits at least the "Iteration"
    // line and the "Converged to energy" line at level 1. Everything
    // above verbosity_ (default 5) should be filtered out.
    logger_solver.set(std::string("verbosity"), 5);
    logger_solver.set(std::string("maximum_iterations"), 3);
    logger_solver.run();

    bool saw_iteration = false, saw_converged = false;
    for (const auto & r : captured) {
      if (r.msg.find("Iteration") != std::string::npos) saw_iteration = true;
      if (r.msg.find("Converged") != std::string::npos) saw_converged = true;
      // Level-gate invariant: no message can come through above the
      // current verbosity_ threshold.
      if (r.level > 5) {
        std::printf("FAIL: logger got level %d message at verbosity 5: %s\n",
                    r.level, r.msg.c_str());
        ++failures;
      }
    }
    REQUIRE(saw_iteration);
    REQUIRE(saw_converged);

    // Silence the solver and confirm nothing further arrives.
    captured.clear();
    logger_solver.set(std::string("verbosity"), 0);
    logger_solver.run();
    for (const auto & r : captured) {
      // Level 0 messages are still legitimate at verbosity 0
      // (unconditional), but level >= 1 must not slip through.
      if (r.level >= 1) {
        std::printf("FAIL: logger got level %d message at verbosity 0: %s\n",
                    r.level, r.msg.c_str());
        ++failures;
      }
    }

    // Clearing the logger restores the stdout default.
    logger_solver.logger(nullptr);
    REQUIRE(!logger_solver.has_logger());
  }

  // ODA simplex projection. Regression test for a HelFEM report where
  // an SCF started from a density projected between two different FEM
  // bases aborted with "Negative natural occupation numbers", the
  // smallest natural occupation being about -5.7e-12 -- some 26000
  // machine epsilons, far above eigensolver roundoff, i.e. a real
  // simplex overshoot rather than noise.
  //
  // The ODA polytope is {lambda >= 0, sum(lambda) <= 1}, and the trial
  // loop only ever scales candidates down, so nothing in the algorithm
  // wants lambda outside it. But the QP solver enforces the simplex
  // only to the accuracy of its constrained linear solve, and an
  // ill-conditioned reduced Hessian (exactly what a cross-basis
  // projection produces) leaves sum(lambda) above one by of order
  // cond * eps. Any overshoot gives the reference density a negative
  // weight, so the mixed density is no longer positive semidefinite.
  {
    using OpenOrbitalOptimizer::HelperRoutines::project_onto_unit_simplex;
    const double eps = std::numeric_limits<double>::epsilon();

    // The reported failure mode: sum marginally above one.
    {
      Vector<double> lam(2);
      lam << 0.5, 0.5 + 3e-12;
      const double before = lam.sum();
      project_onto_unit_simplex<double>(lam);
      REQUIRE(before > 1.0);              // precondition: really overshoots
      REQUIRE(lam.sum() <= 1.0);
      REQUIRE(lam.minCoeff() >= 0.0);
      // Rescaling preserves the direction, so the ratio is untouched.
      REQUIRE(std::abs(lam(0)/lam(1) - 0.5/(0.5 + 3e-12)) < 1e-12);
    }
    // A point strictly inside the simplex must be left alone.
    {
      Vector<double> lam(3);
      lam << 0.25, 0.25, 0.25;
      project_onto_unit_simplex<double>(lam);
      REQUIRE(std::abs(lam(0) - 0.25) <= eps);
      REQUIRE(std::abs(lam(1) - 0.25) <= eps);
      REQUIRE(std::abs(lam(2) - 0.25) <= eps);
    }
    // A point exactly on the sum-cap face is inside, not outside.
    {
      Vector<double> lam(2);
      lam << 0.5, 0.5;
      project_onto_unit_simplex<double>(lam);
      REQUIRE(std::abs(lam.sum() - 1.0) <= eps);
    }
    // Negative entries are clamped, and clamping alone can bring the
    // sum back inside without any rescaling.
    {
      Vector<double> lam(2);
      lam << 1.2, -0.4;
      project_onto_unit_simplex<double>(lam);
      REQUIRE(lam.minCoeff() >= 0.0);
      REQUIRE(lam.sum() <= 1.0 + eps);
    }
    // Gross overshoot rescales rather than truncates: every entry
    // stays positive and the sum lands on the cap.
    {
      Vector<double> lam(3);
      lam << 1.0, 0.5, 0.5;
      project_onto_unit_simplex<double>(lam);
      REQUIRE(std::abs(lam.sum() - 1.0) <= 4*eps);
      REQUIRE(lam.minCoeff() > 0.0);
      REQUIRE(std::abs(lam(0) - 0.5) <= 4*eps);
    }
  }

  std::printf("%s: %d failure(s)\n", __FILE__, failures);
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
