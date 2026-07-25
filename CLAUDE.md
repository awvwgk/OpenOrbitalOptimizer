# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

OpenOrbitalOptimizer is a header-only C++17 library for orbital optimization in quantum chemistry (Hartree-Fock, DFT, and related SCF methods). The library itself is problem-agnostic: it solves the fixed-point equation `FC = CE` in an orthonormal basis, and the caller supplies a Fock builder. Reference paper: J. Phys. Chem. A 129, 5651 (2025).

The library is templated on two types: `SCFSolver<Torb, Tbase>` where `Torb` is the orbital coefficient type (real or complex) and `Tbase` is the (always real) type used for orbital energies and occupations. Valid pairs: `<float,float>`, `<double,double>`, `<std::complex<float>,float>`, `<std::complex<double>,double>`.

## Layout

- `openorbitaloptimizer/scfsolver.hpp` — the entire SCF solver as one header (~4800 lines). All public API lives in `class SCFSolver`; the `public:` section starts around line 3440. Includes DIIS/EDIIS/ADIIS history mixing, the optimal-damping (ODA) polytope step, preconditioned CG / L-BFGS orbital rotations, and Aufbau occupation logic across arbitrary numbers of particle types and symmetry blocks.
- `openorbitaloptimizer/types.hpp` — library typedefs (`Matrix<T>`, `Vector<T>`, `IndexVector`, `FockBuilder`, and the per-block container shorthands) over Eigen.
- `openorbitaloptimizer/eigen_compat.hpp` — small inline helpers filling gaps in Eigen's API (`vectorise_real_imag`, `find_indices_where`, `sort_index_ascending`, `has_nan`/`has_inf`, `expm_antihermitian`).
- `openorbitaloptimizer/quad_support.hpp` — opt-in `_Float128` glue: `numeric_limits`, math overloads and `Eigen::NumTraits`. Only include it if you want quad precision.
- `openorbitaloptimizer/armadillo_compat.hpp` — opt-in compatibility shim exposing the pre-Eigen Armadillo-typed API, kept so downstream Armadillo-only codes (e.g. ERKALE) need not migrate. The core library does not depend on Armadillo.
- `tests/atomtest.cpp` — the main functional test: an atomic SCF/DFT driver using a radial grid (IntegratorXX), Libxc functionals, and BSE-format JSON or ADF-format STO basis sets. Has restricted, unrestricted, and nuclear-electronic-orbital (NEO) drivers. CLI parsed via `tests/cmdline.h`.
- `tests/atomicsolver.hpp` — radial basis abstractions (GTO + STO) used only by `atomtest`.
- `tests/settings_roundtrip.cpp` — runtime coverage of the string-keyed settings façade: `options()` catalog round-trip, unknown/wrong-type key rejection, `print_settings`, citation, and the log-sink callback.
- `tests/{float_real,float_complex,double_complex,quad_real,quad_complex}.cpp` — compile-only template instantiation tests for the non-default `(Torb,Tbase)` pairs.
- `tests/settings_api_check.hpp` — included by every instantiation test; forces instantiation of the settings-façade member templates, which explicit class instantiation does not cover.
- `cmake/` — `OpenOrbitalOptimizerConfig.cmake.in` and an Armadillo-target healing helper for Conda Windows builds.

## Build and test

The library is header-only; "build" really means building the test suite.

```bash
cmake -S . -B objdir -DCMAKE_BUILD_TYPE=Release
cmake --build objdir
ctest --output-on-failure --test-dir objdir
```

A pre-existing `objdir/` is checked into the working tree and is the conventional build directory.

Test dependencies (only required when `OpenOrbitalOptimizer_BUILD_TESTING=ON`, which is the default for top-level builds): Armadillo (always required), Libxc, IntegratorXX, nlohmann_json. CI installs these via conda-forge.

CTest targets defined in `tests/CMakeLists.txt`:
- `openorbopt/atomtest/build` — builds `openorbopt-atomtest`.
- `openorbopt/atomtest/run1` — closed-shell oxygen with PBE/cc-pVDZ.
- `openorbopt/atomtest/run2` — open-shell oxygen (M=3) with PBE/cc-pVDZ.
- `openorbopt/{float-float,cplxfloat-float,cplxdouble-double}/build` — compile-only checks for the alternate template instantiations (these targets are `EXCLUDE_FROM_ALL`).

Run a single test by name, e.g.:
```bash
ctest --test-dir objdir -R openorbopt/atomtest/run2 --output-on-failure
```

Run the atom driver directly (useful when iterating on the solver):
```bash
./objdir/tests/openorbopt-atomtest --Z 8 --M 3 \
  --xfunc GGA_X_PBE --cfunc GGA_C_PBE \
  --basis tests/cc-pvdz.json
```
Other notable flags: `--Q` (charge), `--restricted` (-1=auto), `--Ngrid`, `--sto` (parse ADF STO basis instead of BSE JSON GTO), `--pbasis` (enables NEO mode), `--convthr`, `--lindepthresh`, `--verbosity`.

## Conventions

- Mozilla Public License 2.0; preserve the MPL header on existing files and add it to new ones.
- The solver is a single header. Keep the public API on `SCFSolver` and put helper utilities in the existing namespaces (`OpenOrbitalOptimizer`, `OpenOrbitalOptimizer::HelperRoutines`).
- Solver options are configured through the string-keyed façade — `set(key, value)`, `get_real`/`get_int`/`get_string(key)`, and the static `options()` catalog. There are no per-option typed accessors. When adding a knob, update the member, the matching `set_*`/`get_*` branch, and `options()` together.
- The library does not depend on Libxc/IntegratorXX/nlohmann_json — only the tests do. Do not introduce these (or any other) dependencies into `openorbitaloptimizer/`.
- Use `Tbase(...)` for numeric literals in solver arithmetic; a bare `double` literal silently promotes the computation to double precision in the `float` and `_Float128` instantiations.
- Numeric arguments passed to `log_()` must be cast to `double` explicitly: `_Float128` is not promoted through varargs, so an uncast `Tbase` is undefined behaviour in the quad instantiation. `log_stream_()` has no such restriction.
- The library must remain instantiable for all `(Torb,Tbase)` combinations; when changing template code, build the `openorbopt-instantiation-*` targets to verify. Note that explicit class-template instantiation does not instantiate member templates — `tests/settings_api_check.hpp` exists to cover the settings façade for every pair.
- `objdir/`, `runs/`, `psi4/`, `openorbital.old/`, and various `*~` / `#*#` editor backups in the working tree are local artifacts — do not commit them.
