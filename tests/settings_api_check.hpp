/*
 Copyright (C) 2023- Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

// Compile-time check that the string-keyed settings façade is usable
// from EVERY (Torb, Tbase) instantiation, not just <double, double>.
//
// This exists because the façade previously exposed plain overloads
// `set(const std::string &, Tbase)` and `set(const std::string &, int)`.
// A bare floating literal such as `1e-9` has type `double`, which
// converts to both `int` and a non-double `Tbase` at the same
// (Conversion) rank, so the call was ambiguous and simply did not
// compile for the float and _Float128 instantiations. Likewise
// `set(key, 100u)` was ambiguous for every instantiation. Since the
// typed per-option accessors were removed, `set()` is the only
// configuration path, so that was a hard API break for exactly the
// arbitrary-precision users the library exists to serve.
//
// The class-template instantiation in each test file does NOT cover
// this: `set` is now itself a member template, and explicit class
// instantiation does not instantiate member templates. Only a real
// call does. Hence this header, invoked from every instantiation test.

#include <openorbitaloptimizer/scfsolver.hpp>

#include <complex>
#include <cstddef>
#include <string>

namespace OpenOrbitalOptimizer {
  namespace SettingsApiCheck {

    /// Exercise every argument form the façade must accept. Never
    /// called at run time -- taking its address in main() is enough to
    /// force instantiation and therefore overload resolution.
    template<typename Torb, typename Tbase>
    void check(SCFSolver<Torb, Tbase> & s) {
      // Bare double literal: the case that used to be ambiguous for
      // every Tbase other than double.
      s.set("convergence_threshold", 1e-9);
      // Tbase-typed argument.
      s.set("convergence_threshold", Tbase(1e-9));
      // Unsigned and size_t: used to be ambiguous for all Tbase.
      s.set("maximum_iterations", 100u);
      s.set("maximum_iterations", std::size_t(100));
      // Plain int and bool.
      s.set("maximum_history_length", 8);
      s.set("frozen_occupations", true);
      // String literal and std::string.
      s.set("methods", "DIIS + ODA + CG");
      s.set("error_norm", std::string("rms"));
      // Explicit typed entry points.
      s.set_real("convergence_threshold", Tbase(1e-9));
      s.set_int("verbosity", 0);
      s.set_string("error_norm", "rms");
      // Getters round-trip.
      (void) s.get_real("convergence_threshold");
      (void) s.get_int("verbosity");
      (void) s.get_string("error_norm");
      (void) SCFSolver<Torb, Tbase>::options();
    }

  }
}
