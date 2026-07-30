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

#### New Features
* Convergence is now judged in occupation space as well as in the
  orbital gradient. Two diagnostics measure it, both re-evaluated on
  every read like ``converged``:
    - ``particle_number_error`` -- the largest ``|sum(n) - N|`` over
      the particle types. This is an invariant rather than a
      convergence measure: every step mixes densities that already
      carry the right particle number, so losing any means something
      discarded it, and iterating will not bring it back. It is
      therefore reported but deliberately does not gate convergence,
      since a solver that cannot finish is worse than one that tells
      you its answer is off.
    - ``aufbau_error`` -- the largest occupation sitting above the
      Fermi level or missing from below it, measured against the
      Aufbau filling of the current orbital energies. Orbitals inside
      the Fermi-level degenerate cluster are exempt, that being where
      fractional occupation is legitimate. This one *is* a convergence
      measure and bounds convergence through the new
      ``aufbau_convergence_threshold`` setting (default 1e-6, negative
      to switch the check off).
    - The order matters and is why the bound is a separate predicate,
      ``occupations_converged()``, rather than part of ``converged()``:
      occupation space is put right by the Aufbau cleanup step, which
      needs the gradient to have converged first. A single predicate
      demanding both would block on an error that nothing had yet been
      allowed to fix -- measured on iron, the Aufbau error sits at
      5e-6 and never falls until the cleanup runs. ``run()`` therefore
      orders them gradient, then cleanup, then occupations, and says
      so at the volume of the iteration line if the occupations are
      what is holding the SCF open.
    - Costs nothing on a healthy SCF. On oxygen the Aufbau error is
      already exactly 0 several iterations before the DIIS error
      reaches its threshold, and iron converges in the same 431
      iterations to the same energy with the check on or off.
* New ``"LCIIS"`` method token implementing the least-squares
  commutator in the iterative subspace of Li & Yaron, J. Chem.
  Theory Comput. **12**, 5322 (2016),
  doi:[10.1021/acs.jctc.6b00666](https://doi.org/10.1021/acs.jctc.6b00666).
  Where Pulay's CDIIS minimises ``|| sum_i c_i [F_i, D_i] ||_F^2``,
  LCIIS minimises the commutator between the *predicted* Fock matrix
  and the *predicted* density,
  ``|| [sum_i c_i F_i, sum_j c_j D_j] ||_F^2``, which is genuinely
  quartic in ``c`` and therefore needs the full ``M x M`` grid of
  mixed commutators ``[F_i, D_j]`` rather than just the diagonal.
  The quartic is minimised under ``sum_i c_i = 1`` by Newton's
  method on the Lagrangian, seeded from the CDIIS coefficients.
    - LCIIS is a variant of the extrapolation step, not a step of
      its own: the token implies ``"DIIS"``, so it slots into the
      existing state machine and the A/EDIIS bracketing continues to
      apply. Use it as ``"LCIIS"``, ``"LCIIS + ODA + CG"``, and so on.
      Because the two share the one extrapolation step, asking for
      both (``"DIIS + LCIIS"``) throws rather than silently
      discarding the DIIS request.
    - Any failure of the quartic solve -- singular KKT matrix,
      non-finite iterate, no convergence within the iteration cap, a
      negative target function -- falls back to the CDIIS
      coefficients rather than throwing. LCIIS is a convergence
      accelerator; a bad quartic solve is a reason to take the
      ordinary DIIS step, not to abort the SCF.
    - Three new settings: ``lciis_maximum_history`` (default 6),
      ``lciis_maximum_iterations`` (50) and
      ``lciis_convergence_threshold`` (1e-10). The history cap is
      separate from ``maximum_history_length`` because LCIIS holds
      the whole ``M x M`` commutator grid at once, so its memory
      grows as ``M^2 * n_basis^2`` and its commutator build as
      ``M^2`` rather than DIIS's ``M``.
    - On the oxygen atom with PBE/cc-pVDZ, LCIIS reaches the same
      energies as DIIS in fewer iterations (7 vs 9 at M=1, 8 vs 9 at
      M=3).
* ``atomtest`` gained a ``--methods`` flag that overrides the
  driver's default method mix, so any token combination can be
  exercised from the command line.

#### Breaking Changes
* The default method mix is now ``"DIIS + ODA + LBFGS"`` rather than
  ``"DIIS + ODA + CG"``: the orbital-rotation step is taken by L-BFGS
  instead of preconditioned PR+ CG. Callers that relied on the old
  default can ask for it explicitly with
  ``set("methods", "DIIS + ODA + CG")``. The Armadillo shim's
  ``run()`` default argument and the Python binding's documented
  default follow suit.
* Requesting two methods that occupy the same slot now throws
  instead of silently resolving to one of them. This affects
  ``"CG + LBFGS"``, which previously ran L-BFGS and quietly ignored
  the ``CG`` request -- the two are alternative implementations of
  the single orbital-rotation step, not steps that can both run.
  The new ``"DIIS + LCIIS"`` combination is rejected for the same
  reason. Method strings naming only one of each pair are
  unaffected.

#### Enhancements
* Each setting is now a single self-describing object. An option used
  to be named in three places -- the ``options()`` catalog and one
  branch each in the ``set_*`` and ``get_*`` if-else chains -- so 31
  settings carried 88 copies of their key strings, and adding a knob
  meant editing three lists that nothing checked against each other.
  A catalog entry could disagree with what the dispatch actually
  accepted, silently.
    - A setting's value, key, documentation and validator now live
      together in one declaration: ``Setting<T>`` holds the value and
      derives from an abstract ``SettingBase`` carrying everything
      that does not depend on ``T``. ``options()``, ``set_*`` and
      ``get_*`` all read off the settings themselves, so a key string
      appears exactly once and the catalog cannot drift out of step
      with the dispatch.
    - Settings register themselves as they are constructed, so adding a
      knob is one declaration and nothing else -- there is no list to
      remember to append to, and a setting that is declared but not
      published is no longer possible to write. The registry records
      each setting's offset within the solver rather than its address,
      so the facade keeps working after the solver is moved.
    - ``Setting<T>`` converts implicitly to ``const T &`` and assigns
      from ``T``, so the ~195 uses of these members in the solver
      read as they did before.
    - A setting may carry a write validator or a read recompute, as a
      pointer to the member function implementing it. This covers the
      two settings that validate (``error_norm`` probes the norm,
      ``methods`` parses and canonicalises) and the five read-only
      diagnostics.
    - The typed getters and setters are one shared template body
      differing only in ``T``, replacing three near-identical
      if-else chains.
    - Setting a read-only diagnostic now reports it as read-only
      rather than as an unknown key.
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
* The Aufbau cleanup is no longer rejected on atoms with a fractionally
  occupied Fermi level. It swapped the occupation vector but judged it
  on orbitals that were relaxed for the *old* occupations, so it lost
  to the converged mixed density every time and was discarded --
  silently, its return value being dropped. The occupations reported
  stayed the mixed ones, Rydberg tail and all.
    - The step now relaxes the orbitals at the occupations it chose,
      to stationarity, and compares only the relaxed energy. On iron
      (PBE/cc-pVDZ) that alone took the gap from 1.64e-2 Eh to
      9.74e-4.
    - The skeleton set is established once and held fixed across the
      cleanup's passes. Which orbitals are degenerate is a property of
      the solution being refined, and re-deriving it from each relaxed
      iterate destroyed it: the relaxation moves the cluster apart by
      more than ``optimal_damping_degeneracy_threshold``, after which
      the walk no longer recognises it and hands back a simplex of
      dimension zero.
    - What was left was a curvature the model could not see. ODA
      expands the energy in the occupations at *fixed* orbitals, so it
      misses the response term in
      ``H_relaxed = H_ll - H_lk H_kk^-1 H_kl``, and being too stiff it
      stops short: on iron it put the beta 4s/3d split at 0.826/1.174
      where the SCF finds 0.657/1.343 -- both on the same
      one-parameter line, so this was never a missing degree of
      freedom. Alternating the two halves does not recover it either;
      each is exact given the other, so the alternation has a fixed
      point and sits down at it.
    - Where the Fermi level is a single degenerate pair, the cleanup
      now models the coupling. It relaxes at the sampled points and
      fits a cubic through relaxed data rather than forming the Schur
      complement, which would need orbital response equations the
      solver does not have. The endpoint slopes are free: at an
      orbital-stationary point the orbital-response contribution to
      ``dE/dlambda`` carries a factor ``dE/dkappa = 0``, leaving
      Hellmann-Feynman.
    - Iron's minimum now lands at lambda = 0.6596, against the SCF's
      0.6567, and reaches 1.3e-6 Eh *below* the converged mixed
      density, so the cleanup is accepted and the reported occupations
      are Aufbau. Oxygen accepts at M = 1, 3 and 5 across three method
      mixes with the reference energies unchanged.
    - This stays out of the SCF loop deliberately: it cost about 100
      Fock builds on iron against 576 for the whole run, the loop
      already performs the alternation across its iterations, and the
      free-slope argument holds only at a stationary point.
    - Fermi levels spanning more than two orbitals still fall back to
      the alternation and may still be rejected.
* ``atomtest`` accepts functional id 0. Libxc numbers its functionals
  from one, so a non-positive id means "nothing to add here" -- 0 as
  written by someone wanting exchange only, -1 as returned by
  ``xc_functional_get_number`` for a name it does not know. Only -1 was
  honoured, and only in two of the three builders, so ``--xfunc 0``
  threw and NEO without an electron-proton correlation functional could
  not run at all.
* The occupations reported at convergence are now Aufbau: full below
  the Fermi level, exactly zero above it, fractional only inside the
  degenerate cluster at it. What was reported before was the natural
  occupation vector of a *mixed* density, and a mixture of densities
  carrying different orbitals is not idempotent shell by shell, so a
  nominally full shell came out at ``max_occ - epsilon`` and orbitals
  well above the Fermi level carried ``epsilon``.
    - At convergence the SCF takes one more ODA step with the current
      density left out of the polytope. That is the whole fix, because
      the mixing is the whole problem: every skeleton is an Aufbau
      filling of one common set of orbitals, so a combination of
      skeletons alone has those orbitals as its natural orbitals and
      the combined occupation vector as its occupations, exactly. The
      Aufbau structure is inherited rather than imposed, and the
      Fermi-level fractions come from minimising the energy over the
      skeleton simplex rather than from a filling rule.
    - Leaving the reference out is a reparametrisation rather than a
      new code path: one skeleton is promoted to the ``lambda = 0``
      vertex and dropped from the axes, so ``n`` skeletons are
      described by ``n-1`` parameters -- the simplex they span. The
      polytope is still ``{lambda >= 0, sum(lambda) <= 1}``, so the
      QP, the cubic rays and the backoff scaling are untouched.
    - Largest occupation above the Fermi level at convergence, oxygen
      with PBE/cc-pVDZ: 6.6e-12 to 0 (M=3) and 2.0e-12 to 0 (M=1),
      with the deviation of full shells from ``max_occ`` following.
      Reference energies are unchanged, M=3 by a digit in its favour.
* ``interpolate_density`` no longer hands back a density with
  default-constructed blocks for a particle that has no trial
  occupations; it copies the reference across instead.
* ODA no longer loses particle number. The mixed density's natural
  occupations were snapped to exactly zero below
  ``sqrt(eps) * max_occ`` -- 3.0e-8 for an s block, 2.1e-7 for an f
  block -- with nothing to put the discarded occupation back. Because
  the snapped density is fed forward as the next iterate, the loss
  compounded rather than staying an output artifact, and a converged
  SCF could report occupations summing to less than the requested
  particle number by ~1e-6 electrons on a block with a Rydberg tail.
  Tightening ``convergence_threshold`` did not help: the shortfall is
  a fixed truncation, not an unconverged iterate.
    - The mixed density is a convex combination of same-basis
      densities, so it is positive semidefinite and carries the trace
      its ingredients carried. Both are now enforced -- negative
      occupations are clamped, then the block is rescaled to the
      trace its inputs sum to -- rather than left to the accuracy of
      the eigendecomposition. Nothing is discarded silently.
    - The ``sqrt(eps)`` tolerance was justified by the conditioning
      of a density "projected between basis sets and then mixed", but
      this call site mixes densities built in the same basis in that
      very iteration, so it carries ordinary elementwise roundoff and
      not ``cond * eps``. Measured on the oxygen ODA steps, the
      eigendecomposition's own error is 6.2e-15 while the snapping
      discarded up to 8.1e-9 from a single block.
    - The abort guard is unchanged and still fires below
      ``-eps^(1/4)``; it is checked on the raw occupations, before
      the clean-up could hide anything from it.
    - Measured end to end, the converged occupation sum minus the
      requested particle number goes from -2.5e-12 to 0 (O, M=1),
      -6.2e-12 to 0 (O, M=3), and -6.6e-10 to -1.1e-14 (Fe, M=5).
      ``atomtest`` now checks this, so the existing ``run1`` / ``run2``
      ctest cases fail if it regresses.
* L-BFGS no longer stalls when combined with ODA. On open-shell
  oxygen (PBE/cc-pVDZ, M = 3) ``"ODA + LBFGS"`` used to exhaust
  1000 iterations at −74.1255 Eh, 0.85 Eh above the answer that
  ``"ODA + CG"`` reached in 24; it now converges in 26 to the same
  −74.9722875549 Eh. Three defects contributed, each fixed:
    - The curvature pair recorded the search *direction* ``d`` as
      its ``s``, but the step taken is ``exp(t K)``, a displacement
      of ``t*d``. The paired ``y = g_new - g_old`` was therefore
      measured across a different displacement than ``s``
      described. Since the line search routinely accepts ``t``
      orders of magnitude away from 1 -- 1.4e-8 in the stalling
      run -- the stored pairs were inconsistent with each other by
      the same orders of magnitude, and the resulting two-loop
      direction exceeded the preconditioned-SD one in length by up
      to a factor of 5.7e6.
    - Neither the L-BFGS history nor the PR+ CG direction was
      cleared after an accepted extrapolation, so the recorded
      gradient belonged to the pre-extrapolation iterate and the
      next ``y`` spanned the DIIS jump as well as the rotation.
      ODA accepts already cleared it; DIIS accepts now do too.
      This is what left ``"DIIS + LBFGS"`` needing 200 iterations
      at M = 3, now 15.
    - The two-loop direction replaced the preconditioned-SD one on
      the strength of ``d.g < 0`` alone. Descending is not the same
      as being worth taking: the direction it displaced routinely
      descended several times faster per unit length. The
      replacement now has to retain at least half that rate.
* ODA no longer aborts with "Negative natural occupation numbers"
  on an SCF started from a density projected between two
  different basis sets. Two independent problems fed the same
  crash, and both are fixed.
    - The polytope parameters are now projected onto their
      simplex ``{lambda >= 0, sum(lambda) <= 1}`` before the
      mixed density is built. The bound-constrained solve
      enforces that simplex only to the accuracy of its
      constrained linear solve, so an ill-conditioned reduced
      Hessian -- exactly what a cross-basis projection produces
      -- can leave ``sum(lambda)`` above one by of order
      ``cond * eps``. Any overshoot gives the reference density a
      negative weight, so the mixed density is no longer positive
      semidefinite. The trial loop only ever scales candidates
      *down*, so nothing in the algorithm wanted a step outside
      the simplex in the first place.
    - The natural-occupation noise tolerances are now derived
      from powers of the machine epsilon and scaled by the
      block's maximum occupation: values below ``sqrt(eps)`` are
      clamped to zero, and the abort fires below
      ``-eps^(1/4)``. The old absolute ``10 * eps`` / ``100 * eps``
      cutoffs were right for a freshly built density but far too
      tight for a projected-and-mixed one, whose eigendecomposition
      carries error of order ``cond * eps``. The guard now
      triggers only on occupations negative enough to mean a
      genuinely corrupt density, and the thresholds still tighten
      automatically in ``__float128``.
* ``initialize_with_orbitals`` now resets ``last_oda_via_collapsed_``
  and ``last_polytope_dimension_`` along with the rest of the
  per-run state. Left set from an earlier run on the same solver,
  the former made the convergence-time full-polytope check fire on
  a run that never took a collapsed ODA step.
* The degenerate-cluster walk had two implementations with opposite
  boundaries: the ODA skeleton enumeration included an orbital
  exactly one threshold away from the cluster start, the
  active-rotation count excluded it. The latter sizes the post-ODA
  CG burst *for the clusters the former creates*, and its docstring
  already claimed the two definitions were identical. Both now go
  through one ``degenerate_cluster_end_``, on the ODA boundary.
* Corrected the ADIIS linear-term docstring: the loop actually
  computes ``2 * <D_i - D_0 | F_0>`` (the standard ADIIS model of
  Hu & Yang, JCP 132, 054109), not ``2 * <D_i - D_0 | F_i - F_0>``
  as the old comment claimed. Computation itself is unchanged.

#### Misc.
* Gave each duplicated concept in ``scfsolver.hpp`` a single
  definition. No behaviour change; energies are bit-identical on
  both the ODA and the DIIS paths.
    - The effective convergence threshold had three spellings, two
      of them the same quantity written ``Tbase(0.1)`` and
      ``Tbase(1)/Tbase(10)``. Now
      ``effective_convergence_threshold_()`` and
      ``minimum_useful_descent_()``, which also retires the last
      bare fractional literal in the file.
    - The natural-orbital re-diagonalisation had two copies that had
      already drifted. The drift turns out to be *correct* and is
      now documented rather than removed: ODA mixes convexly, so its
      result must be positive semidefinite, while DIIS extrapolates
      affinely with weights that may be negative, so negative
      occupations are a normal outcome there. ``natural_orbitals_``
      shares the mechanics and takes the noise tolerance as an
      argument; ``require_nonnegative_occupations_`` is the opt-in
      guard. Its doc comment also claimed the sign flip yields
      *increasing* occupations — it yields decreasing ones, which is
      both the convention and what the code wants.
    - The active-natural-orbital extraction (count, allocate, fill,
      tolerance ``10 * max_occ * eps``) was open-coded in both the
      rank-k density build and the rank-k commutator build. Now
      ``active_natural_orbitals_``.
    - ``trace_diff`` open-coded the O(N^2) elementwise trace
      ``tr_of_product_`` already provides, and ``materialise_density``
      was an ODA-local lambda doing per block what
      ``build_density_block_`` does; it is now the member
      ``build_density_blocks_``.
    - The 449-line body of the ODA ``attempts`` loop was never
      re-indented when the loop was wrapped around it. Whitespace
      only — ``git diff -w`` on that commit is empty.
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

