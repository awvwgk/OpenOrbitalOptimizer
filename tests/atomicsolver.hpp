#include <openorbitaloptimizer/scfsolver.hpp>
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace OpenOrbitalOptimizer {
  namespace AtomicSolver {

    enum RadialBasisType {
      STOBASIS,
      GTOBASIS
    };

    class RadialBasis {
    protected:
      /// Angular momentum of the function
      int angular_momentum_;

      inline double Vn(double n, double x) const {
        return std::tgamma(n+1)/std::pow(x,n+1);
      }

      inline double Wn(int n, double x) const {
        return (n-1)/x;
      }

      /// The same closed forms evaluated in an arbitrary scalar type.
      ///
      /// These are the radial integrals themselves, not contractions
      /// with a density, so their precision is a systematic offset
      /// rather than noise. They are templated all the same: an
      /// orthonormal basis built from a double overlap matrix is only
      /// orthonormal to double, and that error does enter every
      /// subsequent quantity as though it were noise.
      template<typename T>
      inline T Vn_t(T n, T x) const {
        using std::tgamma; using std::pow;
        return tgamma(n+T(1))/pow(x,n+T(1));
      }

      template<typename T>
      inline T Wn_t(int n, T x) const {
        return T(n-1)/x;
      }

#if 0
      inline double Enk(int n, int k, double x, int increment=1) const {
        double E=0.0; // E^n_0 = 0
        for(int ik=1;ik<=k;ik+=increment)
          E=ik*(1.0+E)/((n-ik+1.0)*x);
        return E;
      }
#else
      inline double Enk(int n, int k, double x) const {
        double E=0.0; // E^n_0 = 0
        for(int ik=1;ik<=k;ik++)
          E=ik*(1.0+E)/((n-ik+1.0)*x);
        return E;
      }

      inline double binomial(double n, double k) const {
        double binomial = (std::tgamma(n+1.0)/std::tgamma(n-k+1.0))/std::tgamma(k+1.0);
        return binomial;
      }

      template<typename T>
      inline T binomial_t(T n, T k) const {
        using std::tgamma;
        return (tgamma(n+T(1))/tgamma(n-k+T(1)))/tgamma(k+T(1));
      }

      template<typename T>
      inline T Enk_t(T n, T k, T x) const {
        using std::pow;
        T E = T(0);
        for(int j=0;j<(int)k;j++)
          E += binomial_t<T>(n,T(j))*pow(x,T(j));
        E /= binomial_t<T>(n,k)*pow(x,k);
        return E;
      }

      template<typename T>
      inline T Enk_int_t(int n, int k, T x) const {
        T E = T(0);
        for(int ik=1;ik<=k;ik++)
          E = T(ik)*(T(1)+E)/(T(n-ik+1)*x);
        return E;
      }

      inline double Enk(double n, double k, double x) const {
        double E=0.0;
        for(int j=0;j<k;j++)
          E += binomial(n,j)*std::pow(x,j);
        E /= binomial(n,k)*std::pow(x,k);
        return E;
      }
#endif
      /// Type of this radial basis
      RadialBasisType type_;
      /// Coulomb integral
      virtual double Rmnv(int m, int n, int v, double x, double y) const {throw std::logic_error("Not implemented!\n");}

    public:
      RadialBasis(int angular_momentum) : angular_momentum_(angular_momentum) {};
      /// Get radial basis type
      RadialBasisType get_type() const {return type_;}
      /// Evaluate overlap matrix
      virtual Eigen::MatrixXd overlap() const=0;
      /// Evaluate kinetic energy matrix
      virtual Eigen::MatrixXd kinetic(int l) const=0;
      /// Evaluate nuclear attraction matrix
      virtual Eigen::MatrixXd nuclear_attraction() const=0;
      /// Evaluate radial basis functions
      virtual Eigen::MatrixXd eval_f(const Eigen::VectorXd & x) const {throw std::logic_error("Not implemented!\n");};
      /// Evaluate radial basis functions' derivatives
      virtual Eigen::MatrixXd eval_df(const Eigen::VectorXd & x) const {throw std::logic_error("Not implemented!\n");};
      /// Build Coulomb matrix
      virtual Eigen::MatrixXd coulomb(const std::shared_ptr<const RadialBasis> & other, const Eigen::MatrixXd & Pother) const=0;
      /// Return number of basis functions
      virtual size_t nbf() const=0;
    };

    class STOBasis : public RadialBasis {
    private:
      /// STO exponents
      Eigen::VectorXd zeta_;
      /// Tabulated two-electron integrals, per partner shell.
      mutable std::map<const void *, Eigen::MatrixXd> coulomb_tables_;
      /// STO principal quantum numbers
      Eigen::VectorXi n_;

      /// Evaluate two-electron integral
      inline double Rmnv(int m, int n, int v, double x, double y) const {
        double value = std::tgamma(m+n)/(x*y*std::pow(x+y,m+n-1))*(1+Enk(m+n-1,n-v-1,y/x)+Enk(m+n-1,m-v-1,x/y));
        //printf("Rmnv(%i,%i,%i,%e,%e) = %e\n",m,n,v,x,y,value);
        return value;
      }
      /// Pairs of basis functions
      std::vector<std::pair<size_t,size_t>> basis_function_pairs() const {
        std::vector<std::pair<size_t,size_t>> list;
        for(Eigen::Index i=0;i<zeta_.size();i++)
          for(Eigen::Index j=0;j<=i;j++)
            list.push_back(std::make_pair(static_cast<size_t>(i),static_cast<size_t>(j)));
        return list;
      }
    public:
      /// Constructor
      STOBasis(const Eigen::VectorXd & zeta, const Eigen::VectorXi & n, int angular_momentum) : RadialBasis(angular_momentum), zeta_(zeta), n_(n) {
        type_ = STOBASIS;
        if(zeta_.size() != n_.size())
          throw std::logic_error("zeta and n need to be of same size!\n");
      }
      /// Evaluate overlap matrix
      Eigen::MatrixXd overlap() const override {
        auto list(basis_function_pairs());
        Eigen::MatrixXd S(zeta_.size(), zeta_.size());
#pragma omp parallel for
        for(auto pair: list) {
          size_t i=pair.first;
          size_t j=pair.second;
          // Pitzer, page 244
          S(i,j) = S(j,i) = Vn(n_(i)+n_(j), zeta_(i)+zeta_(j));
          //printf("S(%i,%i) %e\n",i,j,S(i,j));
        }
        return S;
      };
      /// Evaluate kinetic energy matrix
      Eigen::MatrixXd kinetic(int am) const override {
        auto list(basis_function_pairs());
        Eigen::MatrixXd T(zeta_.size(), zeta_.size());
#pragma omp parallel for
        for(auto pair: list) {
          size_t i=pair.first;
          size_t j=pair.second;
          // Pitzer, page 244
          T(i,j) = T(j,i) = 0.5*zeta_(i)*zeta_(j)*(
                                                   Wn(n_(i)-am, zeta_(i)) * Wn(n_(j)-am, zeta_(j)) * Vn(n_(i)+n_(j)-2, zeta_(i)+zeta_(j))
                                                   - (Wn(n_(i)-am, zeta_(i)) + Wn(n_(j)-am, zeta_(j))) * Vn(n_(i)+n_(j)-1, zeta_(i)+zeta_(j))
                                                   + Vn(n_(i)+n_(j), zeta_(i) + zeta_(j))
                                                   );
          //printf("T(%i,%i) %e\n",i,j,T(i,j));
        }
        return T;
      }
      /// Evaluate nuclear attraction matrix
      Eigen::MatrixXd nuclear_attraction() const override {
        auto list(basis_function_pairs());
        Eigen::MatrixXd V(zeta_.size(), zeta_.size());
#pragma omp parallel for
        for(auto pair: list) {
          size_t i=pair.first;
          size_t j=pair.second;
          // Pitzer, page 244
          V(i,j) = V(j,i) = Vn(n_(i)+n_(j)-1, zeta_(i)+zeta_(j));
          //printf("V(%i,%i) %e\n",i,j,V(i,j));
        }
        return V;
      }
      /// Evaluate basis functions
      Eigen::MatrixXd eval_f(const Eigen::VectorXd & x) const override {
        Eigen::MatrixXd f = Eigen::MatrixXd::Zero(x.size(), zeta_.size());
#pragma omp parallel for collapse(2)
        for(Eigen::Index iz=0; iz<zeta_.size(); iz++) {
          for(Eigen::Index ix=0; ix<x.size(); ix++) {
            f(ix,iz) = std::pow(x(ix),n_(iz)-1) * std::exp(-zeta_(iz)*x(ix));
          }
        }
        return f;
      }
      /// Evaluate basis functions
      Eigen::MatrixXd eval_df(const Eigen::VectorXd & x) const override {
        Eigen::MatrixXd df = Eigen::MatrixXd::Zero(x.size(), zeta_.size());
#pragma omp parallel for collapse(2)
        for(Eigen::Index iz=0; iz<zeta_.size(); iz++) {
          for(Eigen::Index ix=0; ix<x.size(); ix++) {
            df(ix,iz) = -zeta_(iz) * std::pow(x(ix),n_(iz)-1) * std::exp(-zeta_(iz)*x(ix));
            if(n_(iz)>1) {
              df(ix,iz) += (n_(iz)-1) * std::pow(x(ix),n_(iz)-2) * std::exp(-zeta_(iz)*x(ix));
            }
          }
        }
        return df;
      }
      /// Evaluate Coulomb matrix
      /// Two-electron radial integrals for this pair of shells,
      /// tabulated on first use.
      ///
      /// They depend only on the exponents and quantum numbers, so
      /// they are the same on every Fock build, and evaluating them
      /// inside the contraction meant recomputing them every time --
      /// on a lanthanum atom in AHGBS-7 that is 7.7 million tgamma
      /// calls per build, and a profile put 94 per cent of the whole
      /// calculation in this function with four fifths of that inside
      /// the gamma function. Tabulated, the contraction below is a
      /// matrix-vector product instead.
      ///
      /// Keyed on the address of the other basis, which is sound here
      /// because the bases are built once and outlive every solver
      /// that uses them.
      const Eigen::MatrixXd & coulomb_table(const STOBasis & other) const {
        auto it = coulomb_tables_.find(&other);
        if(it != coulomb_tables_.end()) return it->second;

        auto list(basis_function_pairs());
        const Eigen::Index nother = other.zeta_.size();
        Eigen::MatrixXd R(list.size(), nother*nother);
        for(size_t p=0;p<list.size();p++) {
          const size_t i=list[p].first, j=list[p].second;
          for(Eigen::Index l=0;l<nother;l++)
            for(Eigen::Index k=0;k<nother;k++)
              // Column-major flattening, to match Pother.data().
              R(p, k + l*nother) =
                Rmnv(n_(i)+n_(j), other.n_(k)+other.n_(l), 0,
                     zeta_(i)+zeta_(j), other.zeta_(k)+other.zeta_(l));
        }
        return coulomb_tables_.emplace(&other, std::move(R)).first->second;
      }

      Eigen::MatrixXd coulomb(const STOBasis & other, const Eigen::MatrixXd & Pother) const {
        const Eigen::MatrixXd & R = coulomb_table(other);
        Eigen::Map<const Eigen::VectorXd> p(Pother.data(), Pother.size());
        const Eigen::VectorXd contracted = R * p;

        auto list(basis_function_pairs());
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(zeta_.size(), zeta_.size());
        for(size_t idx=0;idx<list.size();idx++) {
          const size_t i=list[idx].first, j=list[idx].second;
          J(i,j) = J(j,i) = contracted(idx);
        }
        return J;
      }
      /// Radial integrals in an arbitrary scalar type.
      template<typename T>
      inline T Rmnv_t(int m, int n, int v, T x, T y) const {
        using std::tgamma; using std::pow;
        return tgamma(T(m+n))/(x*y*pow(x+y,T(m+n-1)))
               *(T(1)+Enk_int_t<T>(m+n-1,n-v-1,y/x)+Enk_int_t<T>(m+n-1,m-v-1,x/y));
      }

    public:
      template<typename T>
      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> overlap_t() const {
        auto list(basis_function_pairs());
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> S(zeta_.size(), zeta_.size());
        for(auto pair: list) {
          size_t i=pair.first, j=pair.second;
          S(i,j) = S(j,i) = Vn_t<T>(T(n_(i)+n_(j)), T(zeta_(i))+T(zeta_(j)));
        }
        return S;
      }
      template<typename T>
      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> kinetic_t(int am) const {
        auto list(basis_function_pairs());
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> Tm(zeta_.size(), zeta_.size());
        for(auto pair: list) {
          size_t i=pair.first, j=pair.second;
          const T zi = T(zeta_(i)), zj = T(zeta_(j));
          Tm(i,j) = Tm(j,i) = T(0.5)*zi*zj*(
              Wn_t<T>(n_(i)-am, zi) * Wn_t<T>(n_(j)-am, zj) * Vn_t<T>(T(n_(i)+n_(j)-2), zi+zj)
              - (Wn_t<T>(n_(i)-am, zi) + Wn_t<T>(n_(j)-am, zj)) * Vn_t<T>(T(n_(i)+n_(j)-1), zi+zj)
              + Vn_t<T>(T(n_(i)+n_(j)), zi+zj));
        }
        return Tm;
      }
      template<typename T>
      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> nuclear_attraction_t() const {
        auto list(basis_function_pairs());
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> V(zeta_.size(), zeta_.size());
        for(auto pair: list) {
          size_t i=pair.first, j=pair.second;
          V(i,j) = V(j,i) = Vn_t<T>(T(n_(i)+n_(j)-1), T(zeta_(i))+T(zeta_(j)));
        }
        return V;
      }
      template<typename T>
      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>
      eval_f_t(const Eigen::Matrix<T,Eigen::Dynamic,1> & x) const {
        using std::pow; using std::exp;
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> f =
          Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Zero(x.size(), zeta_.size());
        for(Eigen::Index iz=0; iz<zeta_.size(); iz++)
          for(Eigen::Index ix=0; ix<x.size(); ix++)
            f(ix,iz) = pow(x(ix),T(n_(iz)-1)) * exp(-T(zeta_(iz))*x(ix));
        return f;
      }
      template<typename T>
      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>
      eval_df_t(const Eigen::Matrix<T,Eigen::Dynamic,1> & x) const {
        using std::pow; using std::exp;
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> df =
          Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Zero(x.size(), zeta_.size());
        for(Eigen::Index iz=0; iz<zeta_.size(); iz++)
          for(Eigen::Index ix=0; ix<x.size(); ix++) {
            df(ix,iz) = -T(zeta_(iz)) * pow(x(ix),T(n_(iz)-1)) * exp(-T(zeta_(iz))*x(ix));
            if(n_(iz)>1)
              df(ix,iz) += T(n_(iz)-1) * pow(x(ix),T(n_(iz)-2)) * exp(-T(zeta_(iz))*x(ix));
          }
        return df;
      }
      /// Coulomb contraction accumulated in an arbitrary scalar type.
      ///
      /// The radial integrals stay in double: they depend only on the
      /// basis, so their error is a fixed offset that cancels out of
      /// energy differences. The contraction against the density does
      /// not -- it is summed afresh at every occupation the solver
      /// tries, and it is differences between those that the
      /// occupation refinement has to resolve.
      template<typename T>
      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>
      coulomb_scalar(const STOBasis & other,
                     const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & Pother) const {
        auto list(basis_function_pairs());
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> J =
          Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Zero(zeta_.size(), zeta_.size());
        for(auto pair: list) {
          size_t i=pair.first;
          size_t j=pair.second;
          T Jij = T(0);
          for(Eigen::Index l=0;l<Pother.cols();l++)
            for(Eigen::Index k=0;k<Pother.rows();k++)
              Jij += Pother(k,l) * Rmnv_t<T>(n_(i)+n_(j), other.n_(k)+other.n_(l), 0, T(zeta_(i))+T(zeta_(j)), T(other.zeta_(k))+T(other.zeta_(l)));
          J(i,j) = J(j,i) = Jij;
        }
        return J;
      }
      /// Wrapper for the above
      Eigen::MatrixXd coulomb(const std::shared_ptr<const RadialBasis> & other, const Eigen::MatrixXd & Pother) const override {
        assert(other->get_type() == STOBASIS);
        auto pother = std::dynamic_pointer_cast<const STOBasis>(other);
        return coulomb(*pother, Pother);
      }
      /// Return number of basis functions
      size_t nbf() const override {
        return static_cast<size_t>(zeta_.size());
      }
    };

    class GTOBasis : public RadialBasis {
    private:
      /// GTO exponents
      Eigen::VectorXd zeta_;
      /// Tabulated two-electron integrals, per partner shell.
      mutable std::map<const void *, Eigen::MatrixXd> coulomb_tables_;
      /// GTO principal quantum numbers
      Eigen::VectorXi n_;

      /// Evaluate two-electron integral
      inline double Rmnv(int m, int n, int v, double x, double y) const {
        double value = std::tgamma((m+n-1)/2.0)/(x*y*std::pow(x+y,(m+n-3)/2.0))*(1+Enk(0.5*(m+n-3),0.5*(n-v-2),y/x)+Enk(0.5*(m+n-3),0.5*(m-v-2),x/y))/4.0;
        //printf("Rmnv(%i,%i,%i,%e,%e) = %e\n",m,n,v,x,y,value);
        return value;
      }
      /// Pairs of basis functions
      std::vector<std::pair<size_t,size_t>> basis_function_pairs() const {
        std::vector<std::pair<size_t,size_t>> list;
        for(Eigen::Index i=0;i<zeta_.size();i++)
          for(Eigen::Index j=0;j<=i;j++)
            list.push_back(std::make_pair(static_cast<size_t>(i),static_cast<size_t>(j)));
        return list;
      }
    public:
      /// Constructor
      GTOBasis(const Eigen::VectorXd & zeta, const Eigen::VectorXi & n, int angular_momentum) : RadialBasis(angular_momentum), zeta_(zeta), n_(n) {
        type_ = GTOBASIS;
        if(zeta_.size() != n_.size())
          throw std::logic_error("zeta and n need to be of same size!\n");
      };
      /// Evaluate overlap matrix
      Eigen::MatrixXd overlap() const override {
        auto list(basis_function_pairs());
        Eigen::MatrixXd S(zeta_.size(), zeta_.size());
#pragma omp parallel for
        for(auto pair: list) {
          size_t i=pair.first;
          size_t j=pair.second;
          // Pitzer, page 244. Factor one half is missing in paper
          S(i,j) = S(j,i) = 0.5*Vn(0.5*(n_(i)+n_(j)-1),zeta_(i)+zeta_(j));
        }
        return S;
      };
      /// Evaluate kinetic energy matrix
      Eigen::MatrixXd kinetic(int l) const override {
        auto list(basis_function_pairs());
        Eigen::MatrixXd T(zeta_.size(), zeta_.size());
#pragma omp parallel for
        for(auto pair: list) {
          size_t i=pair.first;
          size_t j=pair.second;
          // Pitzer, page 244. Factor one half is missing in paper
          T(i,j) = T(j,i) = zeta_(i)*zeta_(j)*(
                                               Wn(n_(i)-l,2*zeta_(i)) * Wn(n_(j)-l,2*zeta_(j)) * Vn(0.5*(n_(i)+n_(j)-3),zeta_(i)+zeta_(j))
                                               - (Wn(n_(i)-l,2*zeta_(i)) + Wn(n_(j)-l,2*zeta_(j))) * Vn(0.5*(n_(i)+n_(j)-1),zeta_(i)+zeta_(j))
                                               + Vn(0.5*(n_(i)+n_(j)+1),zeta_(i)+zeta_(j))
                                               );
        }
        return T;
      }
      /// Evaluate nuclear attraction matrix
      Eigen::MatrixXd nuclear_attraction() const override {
        auto list(basis_function_pairs());
        Eigen::MatrixXd V(zeta_.size(), zeta_.size());
#pragma omp parallel for
        for(auto pair: list) {
          size_t i=pair.first;
          size_t j=pair.second;
          // Pitzer, page 244. Factor one half is missing in paper
          V(i,j) = V(j,i) = 0.5*Vn(0.5*(n_(i)+n_(j)-2),zeta_(i)+zeta_(j));
        }
        return V;
      }
      /// Evaluate basis functions
      Eigen::MatrixXd eval_f(const Eigen::VectorXd & x) const override {
        Eigen::MatrixXd f = Eigen::MatrixXd::Zero(x.size(), zeta_.size());
#pragma omp parallel for collapse(2)
        for(Eigen::Index iz=0; iz<zeta_.size(); iz++) {
          for(Eigen::Index ix=0; ix<x.size(); ix++) {
            f(ix,iz) = std::pow(x(ix),n_(iz)-1) * std::exp(-zeta_(iz)*x(ix)*x(ix));
          }
        }
        return f;
      }
      /// Evaluate basis functions
      Eigen::MatrixXd eval_df(const Eigen::VectorXd & x) const override {
        Eigen::MatrixXd df = Eigen::MatrixXd::Zero(x.size(), zeta_.size());
#pragma omp parallel for collapse(2)
        for(Eigen::Index iz=0; iz<zeta_.size(); iz++) {
          for(Eigen::Index ix=0; ix<x.size(); ix++) {
            df(ix,iz) = -2.0 * zeta_(iz) * x(ix) * std::pow(x(ix),n_(iz)-1) * std::exp(-zeta_(iz)*x(ix)*x(ix));
            if(n_(iz)>1) {
              df(ix,iz) += (n_(iz)-1) * std::pow(x(ix),n_(iz)-2) * std::exp(-zeta_(iz)*x(ix)*x(ix));
            }
          }
        }
        return df;
      }
      /// Evaluate Coulomb matrix
      /// Two-electron radial integrals for this pair of shells,
      /// tabulated on first use.
      ///
      /// They depend only on the exponents and quantum numbers, so
      /// they are the same on every Fock build, and evaluating them
      /// inside the contraction meant recomputing them every time --
      /// on a lanthanum atom in AHGBS-7 that is 7.7 million tgamma
      /// calls per build, and a profile put 94 per cent of the whole
      /// calculation in this function with four fifths of that inside
      /// the gamma function. Tabulated, the contraction below is a
      /// matrix-vector product instead.
      ///
      /// Keyed on the address of the other basis, which is sound here
      /// because the bases are built once and outlive every solver
      /// that uses them.
      const Eigen::MatrixXd & coulomb_table(const GTOBasis & other) const {
        auto it = coulomb_tables_.find(&other);
        if(it != coulomb_tables_.end()) return it->second;

        auto list(basis_function_pairs());
        const Eigen::Index nother = other.zeta_.size();
        Eigen::MatrixXd R(list.size(), nother*nother);
        for(size_t p=0;p<list.size();p++) {
          const size_t i=list[p].first, j=list[p].second;
          for(Eigen::Index l=0;l<nother;l++)
            for(Eigen::Index k=0;k<nother;k++)
              // Column-major flattening, to match Pother.data().
              R(p, k + l*nother) =
                Rmnv(n_(i)+n_(j), other.n_(k)+other.n_(l), 0,
                     zeta_(i)+zeta_(j), other.zeta_(k)+other.zeta_(l));
        }
        return coulomb_tables_.emplace(&other, std::move(R)).first->second;
      }

      Eigen::MatrixXd coulomb(const GTOBasis & other, const Eigen::MatrixXd & Pother) const {
        const Eigen::MatrixXd & R = coulomb_table(other);
        Eigen::Map<const Eigen::VectorXd> p(Pother.data(), Pother.size());
        const Eigen::VectorXd contracted = R * p;

        auto list(basis_function_pairs());
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(zeta_.size(), zeta_.size());
        for(size_t idx=0;idx<list.size();idx++) {
          const size_t i=list[idx].first, j=list[idx].second;
          J(i,j) = J(j,i) = contracted(idx);
        }
        return J;
      }
      /// Radial integrals in an arbitrary scalar type.
      template<typename T>
      inline T Rmnv_t(int m, int n, int v, T x, T y) const {
        using std::tgamma; using std::pow;
        return tgamma(T(m+n-1)/T(2))/(x*y*pow(x+y,T(m+n-3)/T(2)))
               *(T(1)+Enk_t<T>(T(m+n-3)/T(2),T(n-v-2)/T(2),y/x)
                    +Enk_t<T>(T(m+n-3)/T(2),T(m-v-2)/T(2),x/y))/T(4);
      }

    public:
      template<typename T>
      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> overlap_t() const {
        auto list(basis_function_pairs());
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> S(zeta_.size(), zeta_.size());
        for(auto pair: list) {
          size_t i=pair.first, j=pair.second;
          S(i,j) = S(j,i) = T(0.5)*Vn_t<T>(T(n_(i)+n_(j)-1)/T(2), T(zeta_(i))+T(zeta_(j)));
        }
        return S;
      }
      template<typename T>
      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> kinetic_t(int am) const {
        auto list(basis_function_pairs());
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> Tm(zeta_.size(), zeta_.size());
        for(auto pair: list) {
          size_t i=pair.first, j=pair.second;
          const T zi = T(zeta_(i)), zj = T(zeta_(j));
          Tm(i,j) = Tm(j,i) = zi*zj*(
              Wn_t<T>(n_(i)-am, T(2)*zi) * Wn_t<T>(n_(j)-am, T(2)*zj) * Vn_t<T>(T(n_(i)+n_(j)-3)/T(2), zi+zj)
              - (Wn_t<T>(n_(i)-am, T(2)*zi) + Wn_t<T>(n_(j)-am, T(2)*zj)) * Vn_t<T>(T(n_(i)+n_(j)-1)/T(2), zi+zj)
              + Vn_t<T>(T(n_(i)+n_(j)+1)/T(2), zi+zj));
        }
        return Tm;
      }
      template<typename T>
      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> nuclear_attraction_t() const {
        auto list(basis_function_pairs());
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> V(zeta_.size(), zeta_.size());
        for(auto pair: list) {
          size_t i=pair.first, j=pair.second;
          V(i,j) = V(j,i) = T(0.5)*Vn_t<T>(T(n_(i)+n_(j)-2)/T(2), T(zeta_(i))+T(zeta_(j)));
        }
        return V;
      }
      template<typename T>
      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>
      eval_f_t(const Eigen::Matrix<T,Eigen::Dynamic,1> & x) const {
        using std::pow; using std::exp;
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> f =
          Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Zero(x.size(), zeta_.size());
        for(Eigen::Index iz=0; iz<zeta_.size(); iz++)
          for(Eigen::Index ix=0; ix<x.size(); ix++)
            f(ix,iz) = pow(x(ix),T(n_(iz)-1)) * exp(-T(zeta_(iz))*x(ix)*x(ix));
        return f;
      }
      template<typename T>
      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>
      eval_df_t(const Eigen::Matrix<T,Eigen::Dynamic,1> & x) const {
        using std::pow; using std::exp;
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> df =
          Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Zero(x.size(), zeta_.size());
        for(Eigen::Index iz=0; iz<zeta_.size(); iz++)
          for(Eigen::Index ix=0; ix<x.size(); ix++) {
            df(ix,iz) = -T(2)*T(zeta_(iz)) * x(ix) * pow(x(ix),T(n_(iz)-1)) * exp(-T(zeta_(iz))*x(ix)*x(ix));
            if(n_(iz)>1)
              df(ix,iz) += T(n_(iz)-1) * pow(x(ix),T(n_(iz)-2)) * exp(-T(zeta_(iz))*x(ix)*x(ix));
          }
        return df;
      }
      /// Coulomb contraction accumulated in an arbitrary scalar type.
      ///
      /// The radial integrals stay in double: they depend only on the
      /// basis, so their error is a fixed offset that cancels out of
      /// energy differences. The contraction against the density does
      /// not -- it is summed afresh at every occupation the solver
      /// tries, and it is differences between those that the
      /// occupation refinement has to resolve.
      template<typename T>
      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>
      coulomb_scalar(const GTOBasis & other,
                     const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & Pother) const {
        auto list(basis_function_pairs());
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> J =
          Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Zero(zeta_.size(), zeta_.size());
        for(auto pair: list) {
          size_t i=pair.first;
          size_t j=pair.second;
          T Jij = T(0);
          for(Eigen::Index l=0;l<Pother.cols();l++)
            for(Eigen::Index k=0;k<Pother.rows();k++)
              Jij += Pother(k,l) * Rmnv_t<T>(n_(i)+n_(j), other.n_(k)+other.n_(l), 0, T(zeta_(i))+T(zeta_(j)), T(other.zeta_(k))+T(other.zeta_(l)));
          J(i,j) = J(j,i) = Jij;
        }
        return J;
      }
      /// Wrapper for the above
      Eigen::MatrixXd coulomb(const std::shared_ptr<const RadialBasis> & other, const Eigen::MatrixXd & Pother) const override {
        assert(other->get_type() == GTOBASIS);
        auto pother = std::dynamic_pointer_cast<const GTOBasis>(other);
        return coulomb(*pother, Pother);
      }
      /// Return number of basis functions
      size_t nbf() const override {
        return static_cast<size_t>(zeta_.size());
      }
    };
  }
}
