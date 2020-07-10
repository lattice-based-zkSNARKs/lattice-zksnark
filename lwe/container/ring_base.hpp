#ifndef __RING_BASE__
#define __RING_BASE__

#include "lwe/container/utils.hpp"
#include "lwe/randomness/prg.hpp"
#include <cstdint>
#include <libff/algebra/fields/bigint.hpp>
#include <libsnark/common/libsnark_serialization.hpp>

namespace libsnark {
    template <typename T, T modulus> class Ring {
    public:
        T value;
        static const constexpr T &mod = modulus;

        static LWERandomness::PseudoRandomGenerator *prg;
        static LWERandomness::DiscreteGaussian *dg;

        Ring() : value() {}
        explicit Ring(const T &value) : value(value) {}
        Ring(const Ring &other) : value(other.value) {}

        inline Ring &operator=(const Ring &o) {
            this->value = o.value;
            return *this;
        }

        inline bool operator==(const Ring &other) const {
            return this->value == other.value;
        }

        inline bool operator!=(const Ring &other) const {
            return this->value != other.value;
        }

        inline Ring &operator+=(const Ring &other) {
            this->value += other.value;
            return *this;
        }

        inline Ring &operator-=(const Ring &other) {
            this->value -= other.value;
            return *this;
        }

        inline Ring &operator*=(const Ring &other) {
            this->value *= other.value;
            return *this;
        }

        inline Ring operator+(const Ring &other) const {
            Ring what(*this);
            what += other;
            return what;
        }

        inline Ring operator-(const Ring &other) const {
            Ring what(*this);
            what -= other;
            return what;
        }

        inline Ring operator*(const Ring &other) const {
            Ring what(*this);
            what *= other;
            return what;
        }

        inline Ring squared() const {
            Ring r(*this);
            r *= r;
            return r;
        }

        static Ring zero() { return Ring(0); }
        static Ring one() { return Ring(1); }

        inline Ring &rescale(__uint128_t modulus_scale, uint64_t p_prime) {
            this->value = this->value & (modulus - 1);
            uint64_t mod_res = this->value % p_prime;
            __uint128_t div_interval =
                            modulus | 0 ? modulus / modulus_scale
                            : modulus_scale & 1
                                ? __uint128_t(-1) / modulus_scale
                                : (__uint128_t(1) << 127) / (modulus_scale / 2),
                        res_0 = this->value / div_interval,
                        res_1 = this->value % div_interval,
                        rescale = res_0 + (res_1 > div_interval >> 1 ? 1 : 0),
                        p_prime_round = rescale / p_prime,
                        p_pivot = p_prime_round * p_prime + mod_res;
            uint64_t diff_pivot =
                p_pivot > rescale ? p_pivot - rescale : rescale - p_pivot;

            if (diff_pivot > p_prime / 2) {
                if (p_pivot < rescale && p_pivot + p_prime < modulus_scale)
                    p_pivot += p_prime;
                else if (p_pivot > rescale && p_pivot >= p_prime)
                    p_pivot -= p_prime;
            }
            this->value = p_pivot;
            return *this;
        }

        static Ring random_element() {
            return Ring(prg->bounded((__uint128_t) modulus));
        }

        template <uint64_t LENGTH>
        static void
        random_element_sequence(std::array<Ring<T, modulus>, LENGTH> &_dest) {
            prg->prg_mem_randomize(_dest);
            auto &T_dest = reinterpret_cast<std::array<T, LENGTH> &>(_dest);
            for (uint64_t i = 0; i < LENGTH; i++)
                _dest[i] = Ring<T, modulus>(T_dest[i]);
        }

        static Ring bounded(__uint128_t bound) {
            return Ring(prg->bounded(bound));
        }

        template <uint64_t LENGTH>
        static void
        bounded_sequence(__uint128_t bound,
                         std::array<Ring<T, modulus>, LENGTH> &_dest) {
            for (uint64_t i = 0; i < LENGTH; i++)
                _dest[i] = Ring::bounded(bound);
        }

        static Ring pm_bounded(__uint128_t bound) {
            return Ring(T(0) + prg->pm_bounded(bound));
        }

        template <uint64_t LENGTH>
        static void
        pm_bounded_sequence(__uint128_t bound,
                            std::array<Ring<T, modulus>, LENGTH> &_dest) {
            for (uint64_t i = 0; i < LENGTH; i++)
                _dest[i] = Ring::pm_bounded(bound);
        }

        static Ring discrete_gaussian() { return Ring(T(0) + dg->sample()); }

        template <uint64_t LENGTH>
        static void discrete_gaussian_sequence(
            std::array<Ring<T, modulus>, LENGTH> &_dest) {
            std::array<uint64_t, LENGTH> rnd_src;
            prg->prg_mem_randomize(rnd_src);
            for (uint64_t i = 0; i < LENGTH; i++) {
                auto bucket =
                    --dg->probability_interval.lower_bound(rnd_src[i]);
                _dest[i] = Ring<T, modulus>(T(0) + bucket->second);
            }
        }

        friend std::ostream &operator<<(std::ostream &o,
                                        const Ring<T, modulus> &p) {
            o << (p.value % modulus);
            return o;
        }
        friend std::istream &operator>>(std::istream &i, Ring<T, modulus> &p) {
            i >> p.value;
            return i;
        }
    };

    template <typename T, T modulus>
    LWERandomness::PseudoRandomGenerator *Ring<T, modulus>::prg;

    template <typename T, T modulus>
    LWERandomness::DiscreteGaussian *Ring<T, modulus>::dg;
}

#endif
