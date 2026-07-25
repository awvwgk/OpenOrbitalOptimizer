/*
 *                This Source Code Form is subject to the
 *                terms of the Mozilla Public License, v. 2.0.
 *                If a copy of the MPL was not distributed
 *                with this file, You can obtain one at
 *                http://mozilla.org/MPL/2.0/.
 *
 *           Copyright (c) 2025 Susi Lehtola
 */
#ifndef OPENORBITALOPTIMIZER_EIGEN_COMPAT_HPP
#define OPENORBITALOPTIMIZER_EIGEN_COMPAT_HPP

#include "types.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace OpenOrbitalOptimizer {

  /// Real component type of a (possibly complex) scalar.
  template <class T>
  using RealOf = typename Eigen::NumTraits<T>::Real;

  /// Stack a vector of column vectors into one long column vector. Replaces
  /// arma::join_cols on Cols.
  template <class T>
  Vector<T> join_columns(const std::vector<Vector<T>> & parts) {
    Index total = 0;
    for (const auto & p : parts) total += p.size();
    Vector<T> out(total);
    Index off = 0;
    for (const auto & p : parts) {
      out.segment(off, p.size()) = p;
      off += p.size();
    }
    return out;
  }

  /// Vectorise a real-valued matrix to a column vector (column-major), with no
  /// real/imag splitting.
  template <class T>
  std::enable_if_t<!Eigen::NumTraits<T>::IsComplex, Vector<T>>
  vectorise_real_imag(const Matrix<T> & M) {
    return Eigen::Map<const Vector<T>>(M.data(), M.size());
  }

  /// Vectorise a complex-valued matrix into a real column vector by stacking
  /// the real part on top of the imaginary part. Mirrors the layout the SCF
  /// solver relies on for real-valued optimisation over complex orbital
  /// rotations.
  template <class T>
  std::enable_if_t<Eigen::NumTraits<T>::IsComplex, Vector<RealOf<T>>>
  vectorise_real_imag(const Matrix<T> & M) {
    using R = RealOf<T>;
    Vector<R> out(2 * M.size());
    // M is column-major, so the data pointer streams down columns; this
    // matches arma's storage layout.
    auto realview = Eigen::Map<const Matrix<R>>(reinterpret_cast<const R*>(M.data()),
                                                2, M.size());
    // realview row 0 is the real parts in memory order, row 1 the imag parts.
    out.head(M.size()) = realview.row(0).transpose();
    out.tail(M.size()) = realview.row(1).transpose();
    return out;
  }

  /// Find every index i where pred(v[i]) is true. Stand-in for
  /// arma::find(some_predicate).
  template <class Vec, class Pred>
  IndexVector find_indices_where(const Vec & v, Pred pred) {
    std::vector<Index> hits;
    hits.reserve(v.size());
    for (Index i = 0; i < v.size(); ++i)
      if (pred(v[i]))
        hits.push_back(i);
    IndexVector out(hits.size());
    for (size_t k = 0; k < hits.size(); ++k)
      out[k] = hits[k];
    return out;
  }

  /// Return the indices that sort v in ascending order (stable). Stand-in for
  /// arma::sort_index.
  template <class T>
  IndexVector sort_index_ascending(const Vector<T> & v) {
    IndexVector idx(v.size());
    std::iota(idx.data(), idx.data() + idx.size(), Index{0});
    std::stable_sort(idx.data(), idx.data() + idx.size(),
                     [&](Index a, Index b) { return v[a] < v[b]; });
    return idx;
  }

  /// True iff M contains a NaN. Eigen has allFinite() but not has_nan().
  template <class Mat>
  bool has_nan(const Mat & M) {
    return M.array().isNaN().any();
  }

  /// True iff M contains an infinity.
  template <class Mat>
  bool has_inf(const Mat & M) {
    return (M.array().isInf()).any();
  }

  /// exp(K) for an anti-Hermitian K = -K^\dagger. Computed via the Hermitian
  /// eigendecomposition of iK: iK = U diag(w) U^\dagger with real w, so
  /// exp(K) = exp(-i * iK) = U diag(exp(-i w)) U^\dagger. Returns a matrix of
  /// the same scalar type as K (so complex-typed K produces a complex result;
  /// for a real K the caller is expected to have used the block trick or wrap
  /// in complex first).
  template <class T>
  Matrix<T> expm_antihermitian(const Matrix<T> & K) {
    using R = RealOf<T>;
    // Build iK. For complex T this is a rotation; for real T we need to widen.
    if constexpr (Eigen::NumTraits<T>::IsComplex) {
      Matrix<T> iK = T(R{0}, R{1}) * K; // multiply by i
      Eigen::SelfAdjointEigenSolver<Matrix<T>> es(iK);
      const auto & U = es.eigenvectors();
      const auto & w = es.eigenvalues();
      Vector<T> phase(w.size());
      for (Index i = 0; i < w.size(); ++i)
        phase[i] = std::exp(T(R{0}, -w[i]));
      return U * phase.asDiagonal() * U.adjoint();
    } else {
      // Real anti-symmetric K: promote to complex so that eigenvalues are real.
      Matrix<std::complex<R>> iK(K.rows(), K.cols());
      for (Index c = 0; c < K.cols(); ++c)
        for (Index r = 0; r < K.rows(); ++r)
          iK(r, c) = std::complex<R>(R{0}, R{1}) * K(r, c);
      Eigen::SelfAdjointEigenSolver<Matrix<std::complex<R>>> es(iK);
      const auto & U = es.eigenvectors();
      const auto & w = es.eigenvalues();
      Vector<std::complex<R>> phase(w.size());
      for (Index i = 0; i < w.size(); ++i)
        phase[i] = std::exp(std::complex<R>(R{0}, -w[i]));
      Matrix<std::complex<R>> C = U * phase.asDiagonal() * U.adjoint();
      // The result is real to within round-off for real K.
      Matrix<T> out(K.rows(), K.cols());
      for (Index c = 0; c < K.cols(); ++c)
        for (Index r = 0; r < K.rows(); ++r)
          out(r, c) = static_cast<T>(C(r, c).real());
      return out;
    }
  }

} // namespace OpenOrbitalOptimizer

#endif
