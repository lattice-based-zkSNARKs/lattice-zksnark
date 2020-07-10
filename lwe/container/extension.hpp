#ifndef __EXTENSION_HPP__
#define __EXTENSION_HPP__

#include "lwe/container/ring_base.hpp"

namespace libsnark {
    template <typename T> class Extension {
    public:
        T c0, c1;

        static T non_residue;
        static Extension multiplicative_generator;
        static Extension root_of_unity;
        static size_t s;

        Extension() : c0(), c1() {}
        explicit Extension(const long &x, const long &y = 0) : c0(x), c1(y) {}
        explicit Extension(const T &x, const T &y = T::zero()) : c0(x), c1(y) {}
        Extension(const Extension &r2) : c0(r2.c0), c1(r2.c1) {}

        inline Extension &operator=(const Extension &o) {
            this->c0 = o.c0;
            this->c1 = o.c1;
            return *this;
        }

        inline bool operator==(const Extension &other) const {
            return this->c0 == other.c0 && this->c1 == other.c1;
        }

        inline bool operator!=(const Extension &other) const {
            return this->c0 != other.c0 || this->c1 != other.c1;
        }

        inline Extension &operator+=(const Extension &other) {
            this->c0 += other.c0;
            this->c1 += other.c1;
            return *this;
        }

        inline Extension operator+(const Extension &other) const {
            Extension what(*this);
            what += other;
            return what;
        }

        inline Extension &operator-=(const Extension &other) {
            this->c0 -= other.c0;
            this->c1 -= other.c1;
            return *this;
        }

        inline Extension operator-(const Extension &other) const {
            Extension what(*this);
            what -= other;
            return what;
        }

        inline Extension &operator*=(const Extension &other) {
            const T A = other.c0, B = other.c1;
            const T a = this->c0, b = this->c1;
            const T aA = this->c0 * A, bB = this->c1 * B;
            this->c0 = aA + non_residue * bB;
            this->c1 = (a + b) * (A + B) - aA - bB;
            return *this;
        }

        inline Extension operator*(const Extension &other) const {
            Extension what(*this);
            what *= other;
            return what;
        }

        template <typename oT, oT o_mod>
        inline Extension<Ring<oT, o_mod>>
        lift_ring_multiply(const Extension<Ring<oT, o_mod>> &other) const {
            Extension<Ring<oT, o_mod>> res;
            const Ring<oT, o_mod> A = other.c0, B = other.c1;
            const T a = this->c0, b = this->c1;
            const Ring<oT, o_mod> aA = a.lift_ring_multiply(A),
                                  bB = b.lift_ring_multiply(B);
            res.c0 = aA + non_residue.lift_ring_multiply(bB);
            res.c1 = a.lift_ring_multiply(B) + b.lift_ring_multiply(A);
            return res;
        }

        inline Extension &operator^=(const unsigned long &pwr) {
            if (pwr == 0) {
                this->c0 = T(1);
                this->c1 = T(0);
                return *this;
            }
            unsigned long _pwr = pwr;
            Extension _res(1, 0);
            while (_pwr > 0) {
                if (_pwr & 1u)
                    _res *= *this;
                _pwr >>= 1u;
                *this *= Extension(*this);
            }
            *this = _res;
            return *this;
        }

        inline Extension &operator^=(const libff::bigint<1> &pwr) {
            *this ^= pwr.as_ulong();
            return *this;
        }

        inline Extension operator^(const unsigned long &pwr) const {
            Extension what(*this);
            what ^= pwr;
            return what;
        }

        inline Extension operator^(const libff::bigint<1> &pwr) const {
            Extension what(*this);
            what ^= pwr.as_ulong();
            return what;
        }

        inline Extension squared() const {
            Extension f(*this);
            f *= f;
            return f;
        }

        inline Extension &invert() {
            const T a = this->c0, b = this->c1;
            const T t0 = a.squared();
            const T t1 = b.squared();
            const T t2 = t0 - non_residue * t1;
            const T t3 = t2.inverse();
            this->c0 = a * t3;
            this->c1 = -(b * t3);
            return *this;
        }

        inline Extension inverse() const {
            Extension _e(*this);
            return (_e.invert());
        }

        inline Extension operator-() const {
            Extension r2(*this);
            r2.c0 = -r2.c0;
            r2.c1 = -r2.c1;
            return r2;
        }

        static Extension zero() { return Extension(0, 0); }
        static Extension one() { return Extension(1, 0); }
        static Extension geometric_generator() {
            return Extension::multiplicative_generator;
        }
        static Extension arithmetic_generator() {
            throw std::runtime_error("Fp 2 has no arithmetic generator");
        }

        template <typename nT, nT o_mod>
        inline Extension &
        project_from(const Extension<Ring<nT, o_mod>> &other, nT md_q = o_mod) {
            this->c0.project_from(other.c0, md_q);
            this->c1.project_from(other.c1, md_q);
            return *this;
        }

        template <typename nT, nT o_mod>
        inline void lift_to(Extension<Ring<nT, o_mod>> &other) const {
            this->c0.lift_to(other.c0);
            this->c1.lift_to(other.c1);
        }

        inline Extension &rescale(__uint128_t modulus_scale, uint64_t p_prime) {
            this->c0.rescale(modulus_scale, p_prime);
            this->c1.rescale(modulus_scale, p_prime);
            return *this;
        }

        static Extension random_element() {
            Extension what;
            what.c0 = T::random_element();
            what.c1 = T::random_element();
            return what;
        }

        template <uint64_t LENGTH>
        static void
        random_element_sequence(std::array<Extension<T>, LENGTH> &_dest) {
            auto &flat_T = reinterpret_cast<std::array<T, LENGTH * 2> &>(_dest);
            T::random_element_sequence(flat_T);
        }

        static Extension bounded(__uint128_t bound) {
            Extension what;
            what.c0 = T::bounded(bound);
            what.c1 = T::bounded(bound);
            return what;
        }

        template <uint64_t LENGTH>
        static void bounded_sequence(__uint128_t bound,
                                     std::array<Extension<T>, LENGTH> &_dest) {
            auto &flat_T = reinterpret_cast<std::array<T, LENGTH * 2> &>(_dest);
            T::bounded_sequence(bound, flat_T);
        }

        static Extension pm_bounded(__uint128_t bound) {
            Extension what;
            what.c0 = T::pm_bounded(bound);
            what.c1 = T::pm_bounded(bound);
            return what;
        }

        template <uint64_t LENGTH>
        static void
        pm_bounded_sequence(__uint128_t bound,
                            std::array<Extension<T>, LENGTH> &_dest) {
            auto &flat_T = reinterpret_cast<std::array<T, LENGTH * 2> &>(_dest);
            T::pm_bounded_sequence(bound, flat_T);
        }

        static Extension discrete_gaussian() {
            Extension what;
            what.c0 = T::discrete_gaussian();
            what.c1 = T::discrete_gaussian();
            return what;
        }

        template <uint64_t LENGTH>
        static void
        discrete_gaussian_sequence(std::array<Extension<T>, LENGTH> &_dest) {
            auto &flat_T = reinterpret_cast<std::array<T, LENGTH * 2> &>(_dest);
            T::discrete_gaussian_sequence(flat_T);
        }

        friend std::ostream &operator<<(std::ostream &o,
                                        const Extension<T> &p) {
            o << p.c0 << " " << p.c1;
            return o;
        }
        friend std::istream &operator>>(std::istream &i, Extension<T> &p) {
            i >> p.c0 >> p.c1;
            return i;
        }
    };

    template <typename T> T Extension<T>::non_residue;

    template <typename T> Extension<T> Extension<T>::multiplicative_generator;

    template <typename T> Extension<T> Extension<T>::root_of_unity;

    template <typename T> size_t Extension<T>::s;
}

#endif
