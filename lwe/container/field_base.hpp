#ifndef __FIELD_BASE__
#define __FIELD_BASE__

#include "ring_base.hpp"

namespace libsnark {
    template <typename T, T modulus> class Field {
    public:
        T value;
        static Field multiplicative_generator;
        static Field root_of_unity;
        static size_t s;
        static LWERandomness::PseudoRandomGenerator *prg;
        static LWERandomness::DiscreteGaussian *dg;

        Field() : value{} {}
        explicit Field(const long &x)
            : value(x % modulus + ((x < 0) ? modulus : 0)) {}
        Field(const Field &o) : value(o.value) {}

        inline Field &operator=(const Field &o) {
            this->value = o.value % modulus + ((o.value < 0) ? modulus : 0);
            return *this;
        }

        inline Field &operator=(const T &o) {
            this->value = o % modulus + ((o < 0) ? modulus : 0);
            return *this;
        }

        inline bool operator==(const Field &other) const {
            return this->value == other.value;
        }

        inline bool operator!=(const Field &other) const {
            return this->value != other.value;
        }

        inline Field &operator+=(const Field &other) {
            this->value += other.value;
            this->value %= modulus;
            return *this;
        }

        inline Field operator+(const Field &other) const {
            Field what(*this);
            what += other;
            return what;
        }

        inline Field &operator-=(const Field &other) {
            if (this->value >= other.value)
                this->value -= other.value;
            else
                this->value = modulus - other.value + this->value;
            return *this;
        }

        inline Field operator-(const Field &other) const {
            Field what(*this);
            what -= other;
            return what;
        }

        inline Field &operator*=(const Field &other) {
            this->value *= other.value;
            this->value %= modulus;
            return *this;
        }

        inline Field operator*(const Field &other) const {
            Field what(*this);
            what *= other;
            return what;
        }

        inline Field &operator^=(const unsigned long &pwr) {
            if (pwr == 0) {
                this->value = T(1);
                return *this;
            }
            unsigned long _pwr = pwr;
            Field _res(1);
            while (_pwr > 0) {
                if (_pwr & 1u)
                    _res *= *this;
                _pwr >>= 1u;
                *this *= Field(*this);
            }
            this->value = _res.value;
            return *this;
        }

        inline Field &operator^=(const Field &other) {
            *this ^= ((unsigned long) other.value);
            return *this;
        }

        inline Field &operator^=(const libff::bigint<1> &pwr) {
            *this ^= pwr.as_ulong();
            return *this;
        }

        inline Field operator^(const Field &pwr) const {
            Field what(*this);
            what ^= pwr;
            return what;
        }

        inline Field operator^(const libff::bigint<1> &pwr) const {
            Field what(*this);
            what ^= pwr;
            return what;
        }

        inline Field operator^(unsigned long pwr) const {
            Field what(*this);
            what ^= pwr;
            return what;
        }

        inline Field squared() const {
            Field f(*this);
            f *= f;
            return f;
        }

        inline Field &invert() {
            this->value = LWE::modular_inverse(this->value, modulus);
            return *this;
        }

        inline Field inverse() const {
            Field f(*this);
            return (f.invert());
        }

        inline Field operator-() const {
            Field what(*this);
            what.value = modulus - what.value;
            return what;
        }

        static Field zero() { return Field(0); }
        static Field one() { return Field(1); }
        static Field geometric_generator() {
            return Field::multiplicative_generator;
        }
        static Field arithmetic_generator() { return Field::one(); }

        template <typename nT, nT o_mod>
        Field &project_from(const Ring<nT, o_mod> &o, nT modulus_q = o_mod) {
            __int128_t o_v;
            if (modulus_q & (modulus_q - 1)) {
                o_v = o.value > o_mod >> 1 ? (o.value & (o_mod - 1)) - o_mod
                                           : o.value;
                o_v %= __int128_t(modulus_q);
                if (o_v < 0)
                    o_v += modulus_q;
            } else
                o_v = __int128_t(o.value & (modulus_q - 1));
            if (o_v > __int128_t(modulus_q >> 1))
                o_v -= modulus_q;
            o_v %= __int128_t(modulus);
            if (o_v < 0)
                o_v += modulus;
            this->value = o_v;
            return *this;
        }

        template <typename nT, nT o_mod>
        void lift_to(Ring<nT, o_mod> &other) const {
            other.value = nT(this->value);
        }

        template <typename nT, nT o_mod>
        Ring<nT, o_mod> lift_ring_multiply(const Ring<nT, o_mod> &other) const {
            Ring<nT, o_mod> res(other);
            res.value *= this->value;
            return res;
        }

        static Field random_element() {
            return Field(prg->bounded((__uint128_t) modulus));
        }

        template <uint64_t LENGTH>
        static void
        random_element_sequence(std::array<Field<T, modulus>, LENGTH> &_dest) {
            prg->prg_mem_randomize(_dest);
            auto &T_dest = reinterpret_cast<std::array<T, LENGTH> &>(_dest);
            for (uint64_t i = 0; i < LENGTH; i++)
                _dest[i] = Field<T, modulus>(T_dest[i]);
        }

        friend std::ostream &operator<<(std::ostream &o,
                                        const Field<T, modulus> &p) {
            o << (p.value % modulus);
            return o;
        }
        friend std::istream &operator>>(std::istream &i, Field<T, modulus> &p) {
            i >> p.value;
            return i;
        }
    };

    template <typename T, T modulus>
    Field<T, modulus> Field<T, modulus>::multiplicative_generator;

    template <typename T, T modulus>
    Field<T, modulus> Field<T, modulus>::root_of_unity;

    template <typename T, T modulus> size_t Field<T, modulus>::s;

    template <typename T, T modulus>
    LWERandomness::PseudoRandomGenerator *Field<T, modulus>::prg;

    template <typename T, T modulus>
    LWERandomness::DiscreteGaussian *Field<T, modulus>::dg;
}

#endif
