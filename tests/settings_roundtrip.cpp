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
#include <memory>
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

  // The documented default method mix. Nothing else pins it: the
  // atomtest runs pass --methods, and every check below sets the value
  // before reading it, so a default that drifted away from what the
  // docstring and changelog promise would go unnoticed.
  {
    auto s = make_solver();
    REQUIRE(s.get_string("methods") == "DIIS + ODA + LBFGS");
  }

  // Occupation-space diagnostics. This solver has one block of two
  // orbitals holding two particles, so Aufbau fills the lower orbital
  // and empties the upper one; a converged run has to reproduce that
  // exactly, and to carry the two particles it was asked for.
  {
    auto s = make_solver();
    FockMatrix<double> g(1);
    g[0] = Matrix<double>::Identity(2, 2) * -1.0;
    s.initialize_with_fock(g);
    s.run();
    REQUIRE(s.converged());
    REQUIRE(s.particle_number_error() <= 1e-12);
    REQUIRE(s.aufbau_error() <= 1e-12);
    // Both are also readable through the facade, and re-measure rather
    // than returning a value stored at some earlier point.
    REQUIRE(s.get_real("particle_number_error") == s.particle_number_error());
    REQUIRE(s.get_real("aufbau_error") == s.aufbau_error());
  }

  // A negative Aufbau threshold switches the occupation-space half of
  // the convergence test off, leaving the gradient criterion alone.
  {
    auto s = make_solver();
    FockMatrix<double> g(1);
    g[0] = Matrix<double>::Identity(2, 2) * -1.0;
    s.initialize_with_fock(g);
    REQUIRE(s.get_real("aufbau_convergence_threshold") > 0.0);
    s.set(std::string("aufbau_convergence_threshold"), -1.0);
    REQUIRE(s.occupations_converged());
  }

  // LCIIS method token. LCIIS is a variant of the extrapolation step
  // rather than a step of its own, so the token has to imply DIIS --
  // every state-machine gate gets asked "is DIIS allowed?", never
  // "is LCIIS allowed?".
  {
    auto s = make_solver();
    s.set(std::string("methods"), std::string("lciis"));
    REQUIRE(s.get_string("methods") == "LCIIS");

    s.set(std::string("methods"), std::string("LCIIS + ODA + CG"));
    REQUIRE(s.get_string("methods") == "LCIIS + ODA + CG");

    // Unknown tokens are still rejected, and the diagnostic advertises
    // LCIIS alongside the others.
    bool threw = false;
    std::string what;
    try { s.set(std::string("methods"), std::string("not_a_method")); }
    catch (const std::logic_error & e) { threw = true; what = e.what(); }
    REQUIRE(threw);
    REQUIRE(what.find("LCIIS") != std::string::npos);

    // Asking for both is refused rather than silently resolved: LCIIS
    // substitutes for the CDIIS solve inside the single extrapolation
    // step, so honouring both is not possible and quietly dropping the
    // DIIS request would be worse than an error.
    for (const char * combo : {"DIIS + LCIIS", "LCIIS + DIIS",
                               "diis + lciis + ODA + CG"}) {
      threw = false; what.clear();
      try { s.set(std::string("methods"), std::string(combo)); }
      catch (const std::logic_error & e) { threw = true; what = e.what(); }
      if (!threw) {
        std::printf("FAIL: '%s' was accepted\n", combo);
        ++failures;
      } else {
        REQUIRE(what.find("DIIS") != std::string::npos);
        REQUIRE(what.find("LCIIS") != std::string::npos);
      }
    }
    // A rejected set() must not have disturbed the stored value.
    REQUIRE(s.get_string("methods") == "LCIIS + ODA + CG");

    // The same rule one level down: CG and LBFGS are alternative
    // implementations of the single orbital-rotation step, so asking
    // for both is refused rather than resolved by a silent preference.
    for (const char * combo : {"CG + LBFGS", "LBFGS + CG",
                               "DIIS + ODA + cg + lbfgs"}) {
      threw = false; what.clear();
      try { s.set(std::string("methods"), std::string(combo)); }
      catch (const std::logic_error & e) { threw = true; what = e.what(); }
      if (!threw) {
        std::printf("FAIL: '%s' was accepted\n", combo);
        ++failures;
      } else {
        REQUIRE(what.find("CG") != std::string::npos);
        REQUIRE(what.find("LBFGS") != std::string::npos);
      }
    }
    // Each on its own is still fine.
    s.set(std::string("methods"), std::string("DIIS + ODA + CG"));
    REQUIRE(s.get_string("methods") == "DIIS + ODA + CG");
    s.set(std::string("methods"), std::string("DIIS + ODA + LBFGS"));
    REQUIRE(s.get_string("methods") == "DIIS + ODA + LBFGS");

    // An LCIIS-driven run reaches the same solution as the plain DIIS
    // one on this (trivial, one-block) problem. The point is that the
    // quartic solve runs at all and hands back usable weights rather
    // than falling over -- the real convergence comparison lives in
    // atomtest.
    FockMatrix<double> g(1);
    g[0] = Matrix<double>::Identity(2, 2) * -1.0;

    auto ref = make_solver();
    ref.initialize_with_fock(g);
    ref.set(std::string("methods"), std::string("DIIS"));
    ref.set(std::string("maximum_iterations"), 20);
    ref.run();

    auto lc = make_solver();
    lc.initialize_with_fock(g);
    lc.set(std::string("methods"), std::string("LCIIS"));
    lc.set(std::string("maximum_iterations"), 20);
    lc.run();

    REQUIRE(lc.converged());
    REQUIRE(std::abs(lc.get_energy() - ref.get_energy()) < 1e-10);
  }

  {
    // The settings register themselves with a registry that records
    // offsets, not addresses, so that the facade keeps working after
    // the solver is moved. Had it stored addresses, every lookup below
    // would reach into the moved-from object.
    // The source is heap-allocated and freed immediately after the
    // move, so that a registry of addresses would be left pointing at
    // released storage rather than at a live object that happens to
    // still hold the same values. Keeping the source alive would let
    // the address-based bug pass this test unnoticed -- and short
    // strings survive a move intact under SSO, so reading "methods"
    // back would not have caught it either.
    auto source = std::make_unique<SCFSolver<double, double>>(make_solver());
    source->set(std::string("convergence_threshold"), 1.25e-9);
    source->set(std::string("methods"), std::string("ODA"));
    source->set(std::string("maximum_iterations"), 77);

    SCFSolver<double, double> moved = std::move(*source);
    source.reset();

    REQUIRE(moved.get_real("convergence_threshold") == 1.25e-9);
    REQUIRE(moved.get_string("methods") == "ODA");
    REQUIRE(moved.get_int("maximum_iterations") == 77);

    // Writing through the moved-to object must land in the moved-to
    // object.
    moved.set(std::string("maximum_iterations"), 33);
    REQUIRE(moved.get_int("maximum_iterations") == 33);

    // Every catalogued setting must still resolve in the moved-to
    // object: its registry has to describe all of them, not just the
    // three read back above.
    for (const auto & opt : catalog) {
      const std::string type = opt.type;
      if (type == "real")        (void) moved.get_real(opt.key);
      else if (type == "int")    (void) moved.get_int(opt.key);
      else if (type == "string") (void) moved.get_string(opt.key);
      else { std::printf("unknown type %s\n", opt.type); ++failures; }
    }

    // And it must still solve: the settings the SCF loop reads are the
    // same objects the facade just wrote to.
    FockMatrix<double> g(1);
    g[0] = Matrix<double>::Identity(2, 2) * -1.0;
    moved.initialize_with_fock(g);
    moved.run();
    REQUIRE(moved.converged());
  }

  std::printf("%s: %d failure(s)\n", __FILE__, failures);
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
