# Changelog

<!--
## vX.Y.0 / YYYY-MM-DD (Unreleased)

#### Breaking Changes

#### New Features

#### Enhancements

#### Bug Fixes

#### Misc.
-->


## v0.5.0 / (Unreleased)

#### Enhancements
* ODA now tries a "collapsed" polytope first: when the full
  skeleton enumeration would exceed ``n_particle_types``, each
  particle's skeleton set is replaced by its arithmetic mean --
  one smudged skeleton per particle, one lambda per particle.
  ``npars`` drops from the pathological "high-symmetry guess"
  case (e.g. 28 dimensions in the trace that motivated it) down
  to a handful, cutting vertex evaluations by an order of
  magnitude. If the collapsed pass fails to descend, the full
  skeleton set retries transparently -- trust-region refits and
  everything else apply only to the fallback path. Extra benefit:
  the smudged density has support in every degenerate orbital, so
  the SCF can't collapse onto a symmetry-broken saddle from the
  first-iteration guess.
* Convergence-time full-polytope check. When ``run()`` sees
  ``converged() == true``, ODA is allowed, and the last ODA step
  used the collapsed skeleton set, it invokes
  ``optimal_damping_step(force_full=true)`` before breaking the
  SCF loop. That confirms the converged iterate is stationary in
  the full skeleton set and not just in the collapsed one;
  sub-noise descents (below
  ``0.1 * max(convergence_threshold, K * noise_floor)``, the same
  tolerance the trust-region refit loop uses) are rejected so an
  already-converged SCF doesn't churn at the arithmetic floor.
  After a genuine descent the state machine picks the same
  ``{OrbitalRotation, DIIS, ODA}`` preference the normal post-ODA
  transition uses, so the next step relaxes at the new
  occupations before revisiting DIIS.
* The per-iteration orbital-occupations block in the verbose SCF
  trace now prints before the ``converged()`` check, so the
  occupations at the final iterate reach the log.
* Speed up the ODA quadratic-model construction and the whole
  DIIS / ADIIS / EDIIS extrapolation stack via a coordinated
  refactor of the ``C * diag(n) * C^dagger`` and
  ``(A * B).trace()`` idioms scattered through ``scfsolver.hpp``.
    - New ``build_density_block_`` is the single density-matrix
      construction path. Extracts the occupied natural orbitals
      first, so the outer product is O(k * n_basis^2) instead of
      the naive O(n_basis^3). Now used by ``get_density_matrix_block``
      (so every DIIS-side caller benefits), the ODA polytope-vertex
      cache, ``density_overlap``, the level-shift block, and
      ``interpolate_density``.
    - New ``tr_of_product_`` computes ``tr(A * B)`` in O(n_basis^2)
      via ``(A * B.transpose()).sum()`` instead of the O(n_basis^3)
      matmul-then-trace form. Replaces every ``(A * B).trace()``
      call in the file: ADIIS linear + quadratic, EDIIS quadratic,
      ``density_overlap``, and the ODA gradient / Hessian
      ``trace_diff``.
    - ``adiis_quadratic_term`` and ``ediis_quadratic_term`` now
      pre-cache the block densities so the doubly-nested history
      loop no longer re-materialises the same density
      ``nhist`` times.
    - ``optimal_damping_step`` pre-caches each axis-vertex density
      once via ``materialise_density``; the trace-diff loop then
      does one O(n_basis^2) elementwise trace per call rather than
      an O(n_basis^3) matmul that started from raw ``(C, n)``.
      Total polytope-setup cost drops from ``O(npars^2 * n_basis^3)``
      to ``O(npars * k * n_basis^2 + npars^2 * n_basis^2)`` -- a
      ~20x-30x speed-up already at ``npars = 28`` and increasing
      with basis size, which unlocks useful ODA behaviour on
      heavy-atom / large-basis systems where the older
      implementation appeared to hang after "Roothaan step in
      dimension N".
* Cache DIIS / ADIIS / EDIIS primitives across SCF iterations. Two
  mutable caches keyed on each history entry's stable iteration
  index:
    - ``diis_commutator_cache_`` holds the AO-basis commutator
      ``FP - PF`` per entry per block. Built via the rank-k route
      ``FP = F * (C_occ * diag(n_occ) * C_occ^dagger)`` in
      O(k * n_basis^2), with ``PF = (FP)^dagger`` skipping the
      second full matmul entirely. Dot products of AO-basis
      commutators are unitary-invariant, so the same cache backs
      both ``diis_error_matrix_element`` (no projection) and the
      per-iteration ``diis_residual`` used by ``converged()``.
    - ``trace_DF_cache_`` holds sum-over-blocks ``tr(D_a * F_b)``
      for every entry pair. Every ADIIS / EDIIS linear or
      quadratic entry is one of at most four such primitives, so
      the whole extrapolation-matrix build reduces to map lookups
      once the primitives are populated.
    - ``diis_matrix_cache_`` holds the final DIIS error-matrix
      element ``B(i, j) = -Re(tr(X_i * X_j))`` summed over blocks,
      keyed by the sorted pair of iteration indices. Once cached,
      each subsequent ``B(i, j)`` access is one map lookup with no
      per-block work at all -- the doubly-nested DIIS matrix loop
      only ever visits the row / column of the newest entry.
    - ``density_diff_cache_`` holds the sum-over-blocks Frobenius
      distance ``||D_i - D_j||`` from ``density_matrix_difference``,
      keyed the same way so the ``cleanup()`` sweep runs at
      steady-state cost.
    - Caches are cleared by ``initialize_with_*`` and
      ``reset_history()``. Stale entries otherwise linger, which
      is cheap: one ``Tbase`` per scalar cache and O(n_basis^2)
      per commutator block.
    Per-iteration cost of the DIIS matrix build drops from
    O(nhist^2 * n_basis^3) naive to O(n_basis^2) steady-state
    (only the new entry's row / column populates).
* ODA trust-region refinement: after each ``optimal_damping_step``
  accepts a trial candidate, re-anchor the polytope quadratic
  model at the accepted iterate using the gradient observed there
  (free from the already-computed Fock matrix), re-solve the QP,
  and accept the refined point if it lowers the energy. Exact for
  Hartree-Fock along any linear ray; for DFT, catches the
  residual non-quadraticity the initial axis-vertex data missed.
  Number of refits is capped by the new ``max_oda_refits``
  setting (default 3, 0 disables). Refits pre-check the model's
  predicted improvement against a fraction of the SCF convergence
  threshold and skip the Fock build entirely when the prediction
  is below noise, so idle refits cost nothing on well-conditioned
  problems.
* ODA candidate selection is now model-first. The polytope
  quadratic-model minimum and the 1D cubic Hermite probes (one per
  axis, one per pair-diagonal edge) are ranked by their
  polynomial-predicted energy, and the trial loop evaluates them in
  that order. Along a linear ray in density space the Hartree-Fock
  energy is exactly quadratic and the polynomial fits are exact for
  HF and closely accurate for DFT near convergence, so the sorted
  first candidate is normally the accepted step -- one Fock build
  per ODA call in place of the worst-case
  O(candidates × backoff scales) sweep.
* The 1D probes emit only interior *minima* of the fitted cubic
  (the second root of the derivative is a maximum and is never a
  useful trial step), and de-duplicate the doubled root returned
  when the cubic degenerates to a quadratic. Axis and edge probes
  now share one fitting routine.

#### New Features (continued)
* Route all library log output through a caller-installable
  callback. `SCFSolver::logger(std::function<void(int level, const
  std::string & msg)>)` registers a sink; the library builds each
  formatted message via `vsnprintf` and either hands it to the sink
  (when set) or writes it to stdout (the pre-callback default). The
  `level` argument is the minimum ``verbosity`` at which the
  message would print, so callers can route different severities to
  different destinations. Python bindings expose it as
  `solver.logger(callback)` / `solver.logger(None)` /
  `solver.has_logger()`.
* Every `printf` call in `openorbitaloptimizer/scfsolver.hpp`
  (60 sites, excluding commented-out debug prints) migrated to
  `log_(level, fmt, ...)`. Every `std::cout << ...` matrix / vector
  dump (12 live sites) migrated to `log_stream_(level) << ...`,
  an ostream-style companion proxy that RAII-flushes into the same
  sink. The internal levels reflect the pre-existing
  `if(verbosity_ >= N)` gates; unconditional diagnostic warnings
  and the print_history dump use level 0.


#### Breaking Changes
* The per-option typed getters and setters on `SCFSolver` have been
  removed. Use the string-keyed façade instead:
  `solver.set(key, value)`, `solver.get_real / get_int /
  get_string(key)`, `SCFSolver::options()`. The Armadillo compatibility
  shim keeps its typed forwarders and now routes them through the
  façade.
* `SCFSolver::run()` no longer takes a method-mix argument; it
  consumes the `methods` setting, which is normalised to canonical
  uppercase on `set()`. Migrate `solver.run("ODA + CG")` to
  `solver.set("methods", "ODA + CG"); solver.run()`.
* Python bindings dropped the per-option typed methods (`verbosity`,
  `convergence_threshold`, `maximum_iterations`, ...). Use
  `solver.set(key, value)` or the attribute-style
  `solver.settings.<key> = value`. `SCFSolver` and `OptionInfo` are
  now importable from the `openorbital` package top level.

#### New Features
* SCF convergence threshold is now clamped to an arithmetic-precision
  floor: the effective threshold is
  `max(convergence_threshold, K * noise_floor)`, where `noise_floor`
  is a per-run estimate of the roundoff floor of the DIIS residual
  `C^dagger [F, P] C` frozen from the initial Fock. `K` defaults to
  10 and is tunable via `noise_safety_factor`. `__float128` runs are
  unaffected because their epsilon is tiny; the clamp mainly rescues
  low-precision runs from spinning below what the arithmetic can
  resolve. Callback-driven convergence
  (`callback_convergence_function`) is untouched.
* Introduce a string-keyed settings façade on `SCFSolver`:
  `set(key, value)`, `get_real/get_int/get_string(key)`, and a
  static `options()` catalog listing every knob and read-only
  diagnostic with its type and one-line description. Downstream
  callers now only need to know this triple to reach any setting,
  making JSON/dict-shaped configuration pipe-throughs trivial.
* Add `openorbital.Settings`: an attribute-style proxy that dispatches
  through the C++ catalog, so
  `solver.settings.convergence_threshold = 1e-9` is equivalent to
  `solver.set("convergence_threshold", 1e-9)` with the same catalog
  validation, `dir()` completion, and read-only diagnostic guards.
* Add `SCFSolver::print_settings()` which dumps every catalog entry
  with its current value. Read-only diagnostics that aren't yet
  computable (e.g. `converged` before the first `initialize_with_*`)
  print as `n/a` instead of throwing. Python bindings expose it as
  `solver.print_settings()` and the same output backs
  `str(solver.settings)`.
* Add `SCFSolver::citation()` and `SCFSolver::print_citation()` so
  downstream drivers can echo the canonical reference (Lehtola &
  Burns, J. Phys. Chem. A **129**, 5651 (2025);
  doi:10.1021/acs.jpca.5c02110) without hardcoding it, and add a
  `CITATION.cff` at the repository root so GitHub's "Cite this
  repository" button and citation-tracking tools resolve it
  automatically.
* The three PySCF drivers grew an `options=None` dict argument on
  `kernel()`; entries are forwarded via `solver.set(key, value)`
  before the run.

#### Enhancements
* `L-BFGS` history depth is now controlled by the shared
  `maximum_history_length` setting; the private
  `lbfgs_history_size_` knob has been folded away.
* `brute_force_search_for_lowest_configuration` now saves and
  restores `verbosity` and `frozen_occupations`, so calling it no
  longer permanently silences and thaws the parent solver.
* Removed the deprecated `run_optimal_damping()` alias on the
  Armadillo compatibility shim.

#### Bug Fixes
* Corrected the ADIIS linear-term docstring: the loop actually
  computes ``2 * <D_i - D_0 | F_0>`` (the standard ADIIS model of
  Hu & Yang, JCP 132, 054109), not ``2 * <D_i - D_0 | F_i - F_0>``
  as the old comment claimed. Computation itself is unchanged.

#### Misc.
* Bare ``double`` literals in ``scfsolver.hpp`` arithmetic sites
  (0.0, 1.0, 2.0, 10.0, 100.0, -1.0; 45 occurrences) wrapped in
  ``Tbase(N)``, so ``float`` and ``__float128`` instantiations no
  longer silently promote through ``double``. Bare 0.5 sites
  rewritten as ``Tbase(1)/Tbase(2)`` (7 occurrences), and
  scientific-notation / bare decimal defaults such as
  ``1e-4`` / ``0.02`` wrapped in ``Tbase(...)`` (13 sites) for the
  same reason.


## v0.4.0 / 2026-07-20

#### Breaking Changes
* The linear-algebra backend switched from Armadillo to Eigen 3.4. The public
  API keeps the `SCFSolver<Torb, Tbase>` template signature via a compat shim,
  but callers that reached into Armadillo types directly must migrate to the
  new `OpenOrbitalOptimizer::Matrix<T>` / `Vector<T>` / `IndexVector` aliases.

#### New Features
* [\#38](https://github.com/SusiLehtola/OpenOrbitalOptimizer/pull/38) Port the
  library and `atomtest` from Armadillo to Eigen 3, unlocking arbitrary-precision
  scalar types. Adds a `_Float128` (libquadmath) instantiation of `SCFSolver`
  via `openorbitaloptimizer/quad_support.hpp`.
* [\#10](https://github.com/SusiLehtola/OpenOrbitalOptimizer/pull/10) Bi-level
  optimal-damping + preconditioned CG state machine with a skeleton-density-matrix
  polytope for degenerate shells, orbital-rotation bursts after ODA, and an
  L-BFGS phase. `SCFSolver::run("DIIS + ODA + CG")` is the new default; the
  method mix is user-selectable via `run("DIIS")`, `run("ODA + CG")`, etc.
* [\#10](https://github.com/SusiLehtola/OpenOrbitalOptimizer/pull/10) Optional
  batched Fock-builder callback (`set_batched_fock_builder`) that receives a
  list of trial densities in one call, letting integral / grid setup amortise
  across the ODA polytope axis-vertex sweep.
* [\#10](https://github.com/SusiLehtola/OpenOrbitalOptimizer/pull/10) New
  `--oda` / `--odadegthresh` / `--maxiter` flags in the `atomtest` driver.

#### Enhancements
* Python bindings expose the new API surface: `run(methods=...)`,
  `set_batched_fock_builder`, `has_batched_fock_builder`,
  `clear_batched_fock_builder`, `optimal_damping_degeneracy_threshold`,
  `orbital_rotation_steps_after_oda`, `last_polytope_dimension`,
  `last_active_rotation_count`, `number_of_fock_evaluations`, `converged`.
* The three PySCF drivers (`atomic`, `molecular`, `diatomic`) register
  batched Fock builders automatically and forward `methods` through
  `kernel(methods=...)`.

#### Misc.
* Armadillo is no longer a dependency of the header-only library. `atomtest`
  and the compat shim still require Eigen 3.4+; Armadillo is only pulled in
  by the compat shim's `<-> arma::Mat` conversion helpers when a caller
  explicitly opts in.


## v0.3.0 / 2026-05-21

#### Breaking Changes

#### New Features
* [\#34](https://github.com/SusiLehtola/OpenOrbitalOptimizer/pull/34) Add Python interface!
  This covers the most common `SCFSolver<double, double>` case; others added upon request.
* [\#34](https://github.com/SusiLehtola/OpenOrbitalOptimizer/pull/34) Add a PySCF integration layer
  with three usage-aware drivers, exposed as the `openorbital` package.
* [\#35](https://github.com/SusiLehtola/OpenOrbitalOptimizer/pull/35) Add option `oda_restart_steps`
  to set the number of steps with no DIIS energy improvement after which to use ODA independently of
  DIIS history length. Previously used `maximum_history_length`/2

#### Enhancements

* [\#35](https://github.com/SusiLehtola/OpenOrbitalOptimizer/pull/35) Explicitly symmetrize matices
  for DIIS error calculation to avoid numerical issues.

#### Bug Fixes
* Fix four correctness bugs in SCF solver and CG optimizer
  - get_energy: bounds check used > instead of >=, allowing out-of-bounds access when
    `ihist == orbital_history_.size()`.
  - matricise(vec, dim): the per-block offset was never advanced, so every block read from offset 0
    of the input vector.
  - steepest_descent line search: the "decrease step" branch used std::max, which actually grew the
    step whenever the parabolic prediction was larger than step/20, stalling the search.
  - CG Polak-Ribière: divided by dot(g_prev, g_prev) without guarding against a zero previous gradient,
    propagating NaN into the search direction. Fall back to steepest descent in that case.

#### Misc.
* [\#30](https://github.com/SusiLehtola/OpenOrbitalOptimizer/pull/30) `max(idx)` has been deprecated
  in Armadillo in favor of `index_max()` so switching to the new syntax.
* [\#36](https://github.com/SusiLehtola/OpenOrbitalOptimizer/pull/36) Use modern FindPython with pybind11 module.
* Added a `CLAUDE.md` file to aid agents.



## v0.2.0 / 2025-08-12

#### Breaking Changes

#### New Features

#### Enhancements
 * [\#26](https://github.com/SusiLehtola/OpenOrbitalOptimizer/pull/26) add callback so caller can
   apply its own convergence criteria.
 * MESA should use Aufbau occupations not MOM
 * Undo minimum error criterion by Garza and Scuseria to avoid penalizing large steps leading to a decrease in energy
 * Increase `pure_ediis_factor`
 * Implement callback functionality

#### Bug Fixes
 * [\#27](https://github.com/SusiLehtola/OpenOrbitalOptimizer/pull/27) fix a `Col::subvec()` error
   with minimal basis sets.
 * [\#27](https://github.com/SusiLehtola/OpenOrbitalOptimizer/pull/27) fix UHF with frozen
   occupations by disabling ODA.

#### Misc.
 * Added user guide to readme
 * Deploy docs site
 * Fix GCC 15 warnings


## v0.1.0 / 2025-03-30

#### New Features
 * [\#20](https://github.com/SusiLehtola/OpenOrbitalOptimizer/pull/20) Intf -- allow OOO and
   IntegratorXX to work on Windows.
 * All start-up functionality making library operational.

