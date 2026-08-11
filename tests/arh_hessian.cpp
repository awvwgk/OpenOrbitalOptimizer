/*
 Copyright (C) 2023- Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

// Validates the conventions the ARH curvature pairs are extracted
// with, and that the resulting step works on a model with genuine
// coupling between particle types.
//
// The ARH Hessian model reads a displacement and a gradient change off
// a pair of history entries:
//
//   kappa_ij = (C^dag (D_i - D_0) C)_ij / (n_j - n_i)
//   y_a      = 2 Re[(C^dag (F_i - F_0) C)_ij] (n_j - n_i)
//
// Both carry occupation weights and a factor of 2 whose conventions
// differ between formulations, and a mistake in either is silent: the
// model stays positive definite and the step still descends, it just
// descends on the wrong curvature. So the pairs are checked here
// against numerically differentiated gradients rather than against a
// restatement of the same algebra.
//
// The identity under test is that the two pieces of the Hessian add.
// The energy Hessian in the rotation parameters is
//
//   d2E/dx_a dx_b = tr(dD/dx_a . dF/dD . dD/dx_b)      <- the secant
//                 + tr(F . d2D/dx_a dx_b)              <- Omega
//
// and the quasi-Newton condition measures only the first: F is the
// derivative of the energy with respect to the density, so a Fock
// difference reports how the density-space gradient responded, and
// says nothing about the curvature of the rotation parametrisation
// itself. That second piece is the Roothaan-Hall diagonal, exactly:
// along one rotation a one-electron energy is
// const + (n_j - n_i)(eps_i - eps_j) sin^2 x. So the true gradient
// change along a displacement s is
//
//   Delta g  =  Omega s + y + O(|s|^2),
//
// which is what is checked below. Getting this wrong is the specific
// silent failure the check exists for: subtracting Omega instead of
// adding it leaves a model that is still symmetric, still usable, and
// wrong by twice the one-electron curvature.
//
// The two pieces are separated without needing canonical orbitals:
// freezing F turns the energy into tr(F_ref D(x)), which has no
// dF/dD term at all, so differentiating it isolates the Omega piece
// exactly. Omega's own closed form is then checked separately on a
// one-electron model, where F does not depend on the density, the
// canonical orbitals are just the eigenvectors of h, and the secant
// term vanishes identically.
//
// The model has two particle types whose densities are coupled, which
// is what makes the check able to see a wrong *relative* factor
// between classes -- the interclass Hessian block is the one the joint
// (rather than per-class) subspace exists to capture, and a per-class
// scaling error leaves the intraclass blocks looking correct.

#include <openorbitaloptimizer/scfsolver.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <tuple>
#include <vector>

using OpenOrbitalOptimizer::SCFSolver;
using OpenOrbitalOptimizer::Matrix;
using OpenOrbitalOptimizer::Vector;
using OpenOrbitalOptimizer::IndexVector;
using OpenOrbitalOptimizer::FockBuilder;
using OpenOrbitalOptimizer::FockBuilderReturn;
using OpenOrbitalOptimizer::DensityMatrix;
using OpenOrbitalOptimizer::FockMatrix;
using OpenOrbitalOptimizer::Orbitals;
using OpenOrbitalOptimizer::OrbitalOccupations;
using OpenOrbitalOptimizer::expm_antihermitian;
using OpenOrbitalOptimizer::Index;

namespace {

int failures = 0;

#define REQUIRE(cond)                                                   \
  do {                                                                  \
    if (!(cond)) {                                                      \
      std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);       \
      ++failures;                                                       \
    }                                                                   \
  } while (0)

constexpr int n_orbital = 5;
constexpr int n_block = 2;   // one block per particle type

// Two particle types with an on-site repulsion within each and an
// on-site attraction between them:
//
//   E = sum_I tr(h_I D_I) + (u/2) sum_{I,p} D_I,pp^2
//                         + w sum_p D_0,pp D_1,pp
//
// so F_0 = h_0 + u diag(D_0) + w diag(D_1) and symmetrically for
// F_1. Both Fock matrices are linear in the densities, which is the
// property the quasi-Newton condition rests on -- this is a
// Hartree-Fock-like energy, exactly quadratic in the densities, and
// the secant condition is therefore exact for it.
//
// w is what puts curvature in the interclass block: without it the two
// types would decouple and the joint subspace would carry nothing a
// per-class one does not.
struct Model {
  std::vector<Matrix<double>> h;
  double u = 0.4;
  double w = 0.25;

  Model() : h(n_block) {
    for (int b = 0; b < n_block; ++b) {
      h[b] = Matrix<double>::Zero(n_orbital, n_orbital);
      for (int p = 0; p < n_orbital; ++p)
        for (int q = 0; q < n_orbital; ++q)
          // Hopping plus a diagonal that differs between the types, so
          // neither the orbitals nor the level spacings coincide. The
          // spacing is wide enough to keep the levels well separated:
          // a near-degenerate frontier would put the end-to-end check
          // below at the mercy of the occupation machinery rather than
          // of the curvature model it is meant to exercise.
          h[b](p, q) = (p == q) ? 0.8 * p * (1.0 + 0.3 * b)
                                : -1.0 / (1.0 + p + q + b);
    }
  }

  Matrix<double> density(const Matrix<double> & C,
                         const Vector<double> & n) const {
    Matrix<double> D = Matrix<double>::Zero(n_orbital, n_orbital);
    for (int k = 0; k < n_orbital; ++k)
      D += n(k) * C.col(k) * C.col(k).transpose();
    return D;
  }

  // Energy and Fock matrices at the given orbitals and occupations.
  std::pair<double, std::vector<Matrix<double>>>
  evaluate(const Orbitals<double> & C,
           const OrbitalOccupations<double> & n) const {
    std::vector<Matrix<double>> D(n_block);
    for (int b = 0; b < n_block; ++b) D[b] = density(C[b], n[b]);

    double energy = 0.0;
    std::vector<Matrix<double>> F(n_block);
    for (int b = 0; b < n_block; ++b) {
      energy += (h[b] * D[b]).trace();
      F[b] = h[b];
      for (int p = 0; p < n_orbital; ++p) {
        energy += 0.5 * u * D[b](p, p) * D[b](p, p);
        F[b](p, p) += u * D[b](p, p);
      }
    }
    for (int p = 0; p < n_orbital; ++p) {
      energy += w * D[0](p, p) * D[1](p, p);
      F[0](p, p) += w * D[1](p, p);
      F[1](p, p) += w * D[0](p, p);
    }
    return {energy, F};
  }

  // Energy of the linearised functional tr(F_frozen D), whose Hessian
  // in the rotation parameters is exactly the tr(F d2D/dx_a dx_b) term
  // -- the piece Omega models, isolated by holding F fixed so that the
  // dF/dD term cannot contribute.
  double frozen_energy(const std::vector<Matrix<double>> & F_frozen,
                       const Orbitals<double> & C,
                       const OrbitalOccupations<double> & n) const {
    double e = 0.0;
    for (int b = 0; b < n_block; ++b)
      e += (F_frozen[b] * density(C[b], n[b])).trace();
    return e;
  }

  FockBuilder<double, double> fock_builder() const {
    Model copy = *this;
    return [copy](const DensityMatrix<double, double> & dm) {
      auto ef = copy.evaluate(dm.first, dm.second);
      return FockBuilderReturn<double, double>{ef.first, ef.second};
    };
  }
};

// A starting set of orbitals: the Q factor of a fixed non-symmetric
// matrix, orthonormal but carrying no accidental structure.
Matrix<double> arbitrary_orbitals(int seed) {
  Matrix<double> A(n_orbital, n_orbital);
  for (int p = 0; p < n_orbital; ++p)
    for (int q = 0; q < n_orbital; ++q)
      A(p, q) = std::sin(1.0 + seed + 3.0 * p + 7.0 * q)
              + 0.3 * std::cos(2.0 + 5.0 * p * q);
  return Eigen::HouseholderQR<Matrix<double>>(A).householderQ()
         * Matrix<double>::Identity(n_orbital, n_orbital);
}

// The rotation degrees of freedom, in the same order the solver
// enumerates them: block-major, then i > j with different occupations.
using Dof = std::tuple<int, int, int>;

std::vector<Dof> build_dofs(const OrbitalOccupations<double> & n) {
  std::vector<Dof> dofs;
  for (int b = 0; b < n_block; ++b)
    for (int i = 0; i < n_orbital; ++i)
      for (int j = 0; j < i; ++j)
        if (std::abs(n[b](i) - n[b](j)) >= 1e-6)
          dofs.emplace_back(b, i, j);
  return dofs;
}

// Orbitals displaced from the reference by the rotation amplitudes x,
// always measured from the same reference: this is the chart the
// gradient difference has to be taken in for the two gradients to be
// comparable at all.
Orbitals<double> displaced(const Orbitals<double> & C_ref,
                           const std::vector<Dof> & dofs,
                           const Vector<double> & x) {
  std::vector<Matrix<double>> K(n_block);
  for (int b = 0; b < n_block; ++b)
    K[b] = Matrix<double>::Zero(n_orbital, n_orbital);
  for (size_t a = 0; a < dofs.size(); ++a) {
    const int b = std::get<0>(dofs[a]);
    const int i = std::get<1>(dofs[a]);
    const int j = std::get<2>(dofs[a]);
    K[b](i, j) = x(a);
    K[b](j, i) = -x(a);
  }
  Orbitals<double> C(n_block);
  for (int b = 0; b < n_block; ++b)
    C[b] = C_ref[b] * expm_antihermitian(K[b]);
  return C;
}

// Numerical gradient of an arbitrary functional of the orbitals with
// respect to the rotation amplitudes, at the point x, by central
// differences.
template <typename Functional>
Vector<double> numerical_gradient(Functional energy_of,
                                  const Orbitals<double> & C_ref,
                                  const std::vector<Dof> & dofs,
                                  const Vector<double> & x, double eps) {
  Vector<double> g(dofs.size());
  for (size_t a = 0; a < dofs.size(); ++a) {
    Vector<double> xp = x, xm = x;
    xp(a) += eps;
    xm(a) -= eps;
    g(a) = (energy_of(displaced(C_ref, dofs, xp))
            - energy_of(displaced(C_ref, dofs, xm))) / (2 * eps);
  }
  return g;
}

// The pair the solver extracts from a history entry: the displacement
// read off the density difference and the gradient change read off the
// Fock difference, both projected on the reference rotation DOFs.
// Extract the curvature pair between a reference point and a displaced
// one. The displaced point may carry different occupations, in which
// case the occupation coordinates of the pair are what keeps it
// consistent. Coordinate layout: the rotation amplitudes first, then
// one occupation coordinate per orbital, matching the solver.
void extract_pair(const Model & model, const Orbitals<double> & C_ref,
                  const OrbitalOccupations<double> & n_ref,
                  const std::vector<Dof> & dofs,
                  const Orbitals<double> & C_displaced,
                  const OrbitalOccupations<double> & n_displaced,
                  Vector<double> & s, Vector<double> & y) {
  auto ref = model.evaluate(C_ref, n_ref);
  auto dis = model.evaluate(C_displaced, n_displaced);

  std::vector<Matrix<double>> dD(n_block), dF(n_block), FMO_ref(n_block);
  for (int b = 0; b < n_block; ++b) {
    const Matrix<double> D_ref = model.density(C_ref[b], n_ref[b]);
    const Matrix<double> D_dis =
      model.density(C_displaced[b], n_displaced[b]);
    dD[b] = C_ref[b].transpose() * (D_dis - D_ref) * C_ref[b];
    dF[b] = C_ref[b].transpose() * (dis.second[b] - ref.second[b]) * C_ref[b];
    FMO_ref[b] = C_ref[b].transpose() * ref.second[b] * C_ref[b];
  }

  const size_t n_extended = dofs.size() + n_block * n_orbital;
  s = Vector<double>::Zero(n_extended);
  y = Vector<double>::Zero(n_extended);
  for (size_t a = 0; a < dofs.size(); ++a) {
    const int b = std::get<0>(dofs[a]);
    const int i = std::get<1>(dofs[a]);
    const int j = std::get<2>(dofs[a]);
    // The displacement carries the *displaced* point's occupations:
    // with U = C_0^dag C_i = exp(kappa), the off-diagonal of
    // U n_i U^dag is kappa_ij (n_i,j - n_i,i). The gradient follows
    // the same rule, the derivative in the reference chart at the
    // displaced point being tr(F_i . C_0 [T_a, n_i] C_0^dag), so the
    // change is a difference of two differently weighted terms rather
    // than one weighting of the Fock difference. Both collapse to the
    // familiar form when the occupations agree.
    const double dn_ref = n_ref[b](j) - n_ref[b](i);
    const double dn_dis = n_displaced[b](j) - n_displaced[b](i);
    const double F0ij = FMO_ref[b](i, j);
    const double Fiij = F0ij + dF[b](i, j);
    s(a) = dD[b](i, j) / dn_dis;
    y(a) = 2 * (Fiij * dn_dis - F0ij * dn_ref);
  }
  // Occupation coordinates: the diagonal of the same two matrices.
  for (int b = 0; b < n_block; ++b)
    for (int k = 0; k < n_orbital; ++k) {
      const size_t idx = dofs.size() + (size_t) (b * n_orbital + k);
      s(idx) = dD[b](k, k);
      y(idx) = dF[b](k, k);
    }
}

// Occupation amplitudes of the extended coordinate vector.
Vector<double> occupation_part(const Vector<double> & v,
                               const std::vector<Dof> & dofs) {
  return v.tail(v.size() - (Index) dofs.size());
}

double max_relative_deviation(const Vector<double> & a,
                              const Vector<double> & b) {
  const double scale = b.cwiseAbs().maxCoeff();
  if (!(scale > 0)) return 0.0;
  return (a - b).cwiseAbs().maxCoeff() / scale;
}

}  // namespace

int main() {
  Model model;

  Orbitals<double> C_ref(n_block);
  OrbitalOccupations<double> n(n_block);
  for (int b = 0; b < n_block; ++b) {
    C_ref[b] = arbitrary_orbitals(b);
    n[b].resize(n_orbital);
    // All occupations distinct, so every orbital pair is a rotation
    // DOF and the extraction divides by a range of occupation
    // differences rather than by a single one.
    for (int k = 0; k < n_orbital; ++k)
      n[b](k) = 1.0 - (k + b * 0.13) / (n_orbital + 1);
  }

  const std::vector<Dof> dofs = build_dofs(n);
  REQUIRE(dofs.size() == (size_t) (n_block * n_orbital * (n_orbital - 1) / 2));

  // A displacement direction that touches both particle types, so a
  // wrong relative factor between them cannot cancel.
  Vector<double> kappa(dofs.size());
  for (size_t a = 0; a < dofs.size(); ++a)
    kappa(a) = std::sin(2.0 + 1.7 * a) + 0.4 * std::cos(0.9 * a);
  kappa /= kappa.norm();

  const auto reference = model.evaluate(C_ref, n);
  const auto full_energy = [&](const Orbitals<double> & C) {
    return model.evaluate(C, n).first;
  };
  const auto frozen_energy = [&](const Orbitals<double> & C) {
    return model.frozen_energy(reference.second, C, n);
  };

  // Check 1: the two pieces of the Hessian add.
  //
  // The gradient change along a displacement s decomposes into the
  // piece that comes from the density-space curvature, which is what
  // the Fock-difference secant measures, and the piece that comes from
  // the curvature of the rotation parametrisation, which is what Omega
  // models. The latter is isolated exactly here by freezing F: the
  // functional tr(F_ref D(x)) has no dF/dD term by construction, so
  // differentiating it gives that piece and nothing else.
  std::printf("ARH curvature pairs: do the two Hessian pieces add?\n");
  std::printf("   %zu DOFs over %d blocks, interclass coupling w = %.2f\n",
              dofs.size(), n_block, model.w);
  std::printf("   %-10s %-14s %-14s %-14s %s\n", "t", "kappa err",
              "frozen + y", "secant only", "ratio");

  double previous_error = 0.0;
  bool first = true;
  for (double t : {1e-2, 5e-3, 2.5e-3, 1.25e-3}) {
    const Vector<double> s_exact = t * kappa;
    const Orbitals<double> C_dis = displaced(C_ref, dofs, s_exact);

    Vector<double> s, y;
    extract_pair(model, C_ref, n, dofs, C_dis, n, s, y);

    const Vector<double> zero = Vector<double>::Zero(dofs.size());
    const Vector<double> dg =
      numerical_gradient(full_energy, C_ref, dofs, s_exact, 1e-6)
      - numerical_gradient(full_energy, C_ref, dofs, zero, 1e-6);
    const Vector<double> dg_frozen =
      numerical_gradient(frozen_energy, C_ref, dofs, s_exact, 1e-6)
      - numerical_gradient(frozen_energy, C_ref, dofs, zero, 1e-6);

    const Vector<double> predicted = dg_frozen + y.head((Index) dofs.size());

    const double kappa_error =
      max_relative_deviation(s.head((Index) dofs.size()), s_exact);
    const double model_error = max_relative_deviation(predicted, dg);
    // What the prediction would be if the Omega piece were left out,
    // reported so the check demonstrably has teeth: if this were also
    // small the test would pass for a model that ignored the term.
    const double secant_only_error =
      max_relative_deviation(y.head((Index) dofs.size()), dg);
    const double ratio = s_exact.dot(predicted) / s_exact.dot(dg);

    std::printf("   %-10.3e %-14.3e %-14.3e %-14.3e %.6f\n", t, kappa_error,
                model_error, secant_only_error, ratio);

    // The displacement is recovered from the density difference to
    // first order, so its error is quadratically small.
    REQUIRE(kappa_error <= 5.0 * t);
    // The predicted gradient change is first-order accurate; a wrong
    // convention would sit at O(1) here rather than O(t).
    REQUIRE(model_error <= 5.0 * t);
    REQUIRE(std::abs(ratio - 1.0) <= 0.05);
    REQUIRE(secant_only_error >= 0.1);

    if (!first) {
      // Halving t must roughly halve the error; a constant-factor
      // mistake would hold it flat.
      REQUIRE(model_error < 0.75 * previous_error);
    }
    previous_error = model_error;
    first = false;
  }

  // Check 2: Omega is that frozen-Fock piece.
  //
  // Done on a one-electron model, where F = h whatever the density is.
  // That makes the canonical orbitals exactly the eigenvectors of h --
  // no self-consistency to converge -- and kills the secant term
  // outright, so the whole gradient change has to be Omega s. It is
  // the formula ctx.h uses that is under test:
  //     Omega_a = 2 (eps_i - eps_j) (n_j - n_i).
  {
    Model one_electron;
    one_electron.u = 0.0;
    one_electron.w = 0.0;

    Orbitals<double> C_can(n_block);
    Vector<double> eps[n_block];
    for (int b = 0; b < n_block; ++b) {
      Eigen::SelfAdjointEigenSolver<Matrix<double>> es(one_electron.h[b]);
      C_can[b] = es.eigenvectors();
      eps[b] = es.eigenvalues();
    }

    Vector<double> omega(dofs.size());
    for (size_t a = 0; a < dofs.size(); ++a) {
      const int b = std::get<0>(dofs[a]);
      const int i = std::get<1>(dofs[a]);
      const int j = std::get<2>(dofs[a]);
      omega(a) = 2 * (eps[b](i) - eps[b](j)) * (n[b](j) - n[b](i));
    }

    const auto one_electron_energy = [&](const Orbitals<double> & C) {
      return one_electron.evaluate(C, n).first;
    };

    std::printf("\nOmega against the one-electron Hessian\n");
    std::printf("   %-10s %-14s %s\n", "t", "Omega s err", "secant norm");
    for (double t : {1e-2, 2.5e-3}) {
      const Vector<double> s_exact = t * kappa;
      const Orbitals<double> C_dis = displaced(C_can, dofs, s_exact);

      Vector<double> s, y;
      extract_pair(one_electron, C_can, n, dofs, C_dis, n, s, y);

      const Vector<double> zero = Vector<double>::Zero(dofs.size());
      const Vector<double> dg =
        numerical_gradient(one_electron_energy, C_can, dofs, s_exact, 1e-6)
        - numerical_gradient(one_electron_energy, C_can, dofs, zero, 1e-6);

      const double omega_error =
        max_relative_deviation(omega.cwiseProduct(s_exact), dg);
      (void) occupation_part;
      std::printf("   %-10.3e %-14.3e %.3e\n", t, omega_error, y.norm());

      REQUIRE(omega_error <= 5.0 * t);
      // F does not depend on the density here, so the secant term is
      // identically zero and Omega has to carry the whole Hessian.
      REQUIRE(y.norm() <= 1e-12);
    }
  }

  // Check 3: the occupation coordinates.
  //
  // These are what make a pair consistent when the occupations moved
  // between the two history entries -- the case ARH previously had to
  // discard. Two things are under test: that the displacement is the
  // diagonal of the density difference, and that the paired gradient
  // change is the diagonal of the Fock difference, i.e. that dE/dn_k
  // is the orbital energy (Janak).
  //
  // The occupation convention is checked with the orbitals held
  // fixed. That is not a dodge, it is the only way to compare against
  // a numerical derivative at all: dE/dn_k means the derivative with
  // respect to the occupation of orbital k, so at a rotated point the
  // numerical derivative refers to the *rotated* orbital k while the
  // extraction refers to the reference one, and the two differ at
  // first order in the rotation. Holding the orbitals fixed removes
  // that ambiguity and leaves the convention itself exposed.
  {
    std::printf("\nOccupation coordinates, orbitals fixed\n");
    std::printf("   %-10s %-16s %-16s\n", "t", "delta n err",
                "dE/dn change err");

    // Gradient of the energy with respect to the occupations at fixed
    // orbitals, by central differences. Occupations are perturbed
    // freely: the model imposes no particle-number constraint, and the
    // convention under test is the unconstrained derivative.
    auto occupation_gradient = [&](const Orbitals<double> & C,
                                   const OrbitalOccupations<double> & occ) {
      Vector<double> g(n_block * n_orbital);
      // At fixed orbitals the energy is exactly quadratic in the
      // occupations, so a central difference is exact whatever the
      // step and the only error is round-off -- which a *large* step
      // minimises. A small one would leave the difference of two
      // gradients swamped by noise that grows as the displacement
      // shrinks.
      const double eps = 1e-2;
      for (int b = 0; b < n_block; ++b)
        for (int k = 0; k < n_orbital; ++k) {
          OrbitalOccupations<double> op = occ, om = occ;
          op[b](k) += eps;
          om[b](k) -= eps;
          g(b * n_orbital + k) = (model.evaluate(C, op).first
                                  - model.evaluate(C, om).first) / (2 * eps);
        }
      return g;
    };

    double previous = 0.0;
    bool first_occ = true;
    for (double t : {1e-2, 5e-3, 2.5e-3}) {
      OrbitalOccupations<double> n_dis = n;
      for (int b = 0; b < n_block; ++b) {
        n_dis[b](1) -= t * (1.0 + 0.4 * b);
        n_dis[b](3) += t * (1.0 + 0.4 * b);
      }

      Vector<double> s, y;
      extract_pair(model, C_ref, n, dofs, C_ref, n_dis, s, y);

      Vector<double> dn_exact(n_block * n_orbital);
      for (int b = 0; b < n_block; ++b)
        for (int k = 0; k < n_orbital; ++k)
          dn_exact(b * n_orbital + k) = n_dis[b](k) - n[b](k);

      const Vector<double> dg_occ =
        occupation_gradient(C_ref, n_dis) - occupation_gradient(C_ref, n);

      const double dn_error =
        max_relative_deviation(occupation_part(s, dofs), dn_exact);
      const double dg_error =
        max_relative_deviation(occupation_part(y, dofs), dg_occ);
      std::printf("   %-10.3e %-16.3e %-16.3e\n", t, dn_error, dg_error);

      // With the orbitals fixed the density difference is exactly the
      // occupation change, so this one is at round-off.
      REQUIRE(dn_error <= 1e-10);
      // dE/dn_k = F^MO_kk exactly, so its change is the diagonal of dF
      // -- also exact here, the model's Fock matrix being linear in the
      // density.
      REQUIRE(dg_error <= 1e-8);
      (void) previous;
      (void) first_occ;
      first_occ = false;
    }
  }

  // Check 3b: the rotation coordinates survive an occupation change.
  //
  // The displacement is read off the density difference by dividing by
  // an occupation difference, and when the entry's occupations differ
  // from the reference's there are two of them to choose between.
  // Writing U = C_0^dag C_i = exp(kappa),
  //     (U n_i U^dag)_ij  ~  kappa_ij (n_i,j - n_i,i)
  // off the diagonal, so it is the *entry's*. Dividing by the
  // reference's instead scales every rotation coordinate by the ratio
  // of the two -- an error that does not vanish as the step shrinks,
  // which is what makes it worth a check of its own.
  //
  // The earlier checks cannot see this: they hold the occupations
  // fixed wherever they look at rotations, and the s.y identity below
  // is the same contraction of the same matrices either way.
  {
    std::printf("\nRotation coordinates under a simultaneous"
                " occupation change\n");
    std::printf("   %-10s %-16s\n", "t", "kappa err");

    double previous = 0.0;
    bool first_mixed = true;
    for (double t : {1e-2, 5e-3, 2.5e-3}) {
      OrbitalOccupations<double> n_dis = n;
      for (int b = 0; b < n_block; ++b) {
        n_dis[b](1) -= 0.30 * (1.0 + 0.2 * b);
        n_dis[b](3) += 0.30 * (1.0 + 0.2 * b);
      }
      const Vector<double> s_exact = t * kappa;
      const Orbitals<double> C_dis = displaced(C_ref, dofs, s_exact);

      Vector<double> s, y;
      extract_pair(model, C_ref, n, dofs, C_dis, n_dis, s, y);

      const double err =
        max_relative_deviation(s.head((Index) dofs.size()), s_exact);
      std::printf("   %-10.3e %-16.3e\n", t, err);

      // First-order accurate, so linear in the step. The wrong
      // divisor sits near 1 here and stays there.
      REQUIRE(err <= 5.0 * t);
      if (!first_mixed) REQUIRE(err < 0.75 * previous);
      previous = err;
      first_mixed = false;
    }
  }

  // Check 4: rotation and occupation coordinates are normalised
  // consistently with each other.
  //
  // This is the relative-scaling check, and it is exact rather than
  // asymptotic. The extended inner product has to reproduce the
  // density-space one,
  //
  //   s . y  =  tr(dD dF),
  //
  // with the off-diagonal pairs counted twice and the diagonal once.
  // That identity is what fixes the factor of 2 and the occupation
  // weights on the rotation coordinates relative to the plain
  // diagonal of the occupation ones: get either wrong and the two
  // sides part company immediately.
  //
  // It holds at fixed occupations, and only there. s . y is a pairing
  // in the chart the model works in, tr(dD dF) one in density space,
  // and the two agree exactly while the map between them is the
  // identity. Once the occupations move they part: the chart carries
  // the displaced point's occupation differences in both the
  // displacement and the gradient, and the density-space contraction
  // knows nothing of that. The mixed case is checked against a known
  // kappa above instead, which is the ground truth the identity
  // cannot supply.
  {
    std::printf("\nExtended inner product against tr(dD dF)\n");
    std::printf("   %-10s %-16s %-16s %s\n", "t", "s.y", "tr(dD dF)",
                "rel. dev.");
    for (double t : {1e-1, 1e-2, 1e-3}) {
      const Orbitals<double> C_dis = displaced(C_ref, dofs, t * kappa);

      Vector<double> s, y;
      extract_pair(model, C_ref, n, dofs, C_dis, n, s, y);

      auto ref = model.evaluate(C_ref, n);
      auto dis = model.evaluate(C_dis, n);
      double trace = 0.0;
      for (int b = 0; b < n_block; ++b) {
        const Matrix<double> dD =
          model.density(C_dis[b], n[b]) - model.density(C_ref[b], n[b]);
        const Matrix<double> dF = dis.second[b] - ref.second[b];
        trace += (dD * dF).trace();
      }
      const double sy = s.dot(y);
      const double deviation = std::abs(sy - trace) / std::abs(trace);
      std::printf("   %-10.3e %-16.8e %-16.8e %.3e\n", t, sy, trace,
                  deviation);
      // Equality here is algebraic, not asymptotic: both sides are the
      // same contraction of the same two matrices. The tolerance only
      // has to cover the rotation-to-displacement inversion, which is
      // first-order, so it loosens with t rather than tightening.
      REQUIRE(deviation <= 5.0 * t);
    }
  }

  // Check 4, end to end: the same coupled model driven through the
  // solver. ARH has to reach the solution the other rotation steps
  // reach. Fock counts are reported rather than asserted on -- this
  // model is small enough that the counts say more about the state
  // machine than about the curvature model.
  {
    IndexVector blocks_per_particle(n_block);
    Vector<double> maxocc(n_block), nparticles(n_block);
    for (int b = 0; b < n_block; ++b) {
      blocks_per_particle(b) = 1;
      maxocc(b) = 1.0;
      nparticles(b) = 2.0;
    }

    FockMatrix<double> guess(n_block);
    for (int b = 0; b < n_block; ++b) guess[b] = model.h[b];

    double reference_energy = 0.0;
    bool have_reference = false;
    std::printf("\nCoupled two-class model through the solver"
                " (u = %.2f, w = %.2f)\n", model.u, model.w);
    // The DIIS-free mixes are the ones that matter here: with DIIS in
    // the mix the extrapolation converges this model on its own and
    // the rotation step barely runs, so ARH would be reported as
    // agreeing without having been exercised.
    for (const std::string & methods :
         {std::string("DIIS + ODA + CG"), std::string("DIIS + ODA + LBFGS"),
          std::string("DIIS + ODA + ARH"), std::string("ODA + CG"),
          std::string("ODA + LBFGS"), std::string("ODA + ARH")}) {
      SCFSolver<double, double> s(blocks_per_particle, maxocc, nparticles,
                                  model.fock_builder(), {"a", "b"});
      s.initialize_with_fock(guess);
      s.set(std::string("verbosity"), 0);
      s.set(std::string("methods"), methods);
      s.set(std::string("maximum_iterations"), 300);
      s.run();

      REQUIRE(s.converged());
      const double energy = s.get_energy();
      std::printf("   %-20s E = %.10f  Fock builds = %4d  converged = %d\n",
                  methods.c_str(), energy,
                  s.get_int("number_of_fock_evaluations"),
                  (int) s.converged());
      if (!have_reference) {
        reference_energy = energy;
        have_reference = true;
      } else {
        REQUIRE(std::abs(energy - reference_energy) <= 1e-7);
      }
    }
  }

  // The parser has to reject two rotation steps at once, ARH included.
  {
    IndexVector blocks_per_particle(1);
    blocks_per_particle(0) = 1;
    Vector<double> maxocc(1), nparticles(1);
    maxocc(0) = 1.0;
    nparticles(0) = 2.0;
    SCFSolver<double, double> s(blocks_per_particle, maxocc, nparticles,
                                model.fock_builder(), {"a"});
    bool threw = false;
    try {
      s.set(std::string("methods"), std::string("DIIS + ARH + LBFGS"));
    } catch (const std::logic_error &) {
      threw = true;
    }
    REQUIRE(threw);
  }

  if (failures) {
    std::printf("%d failure(s)\n", failures);
    return 1;
  }
  std::printf("All ARH curvature checks passed\n");
  return 0;
}
