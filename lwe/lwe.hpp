#ifndef __LWE_HPP__
#define __LWE_HPP__

#include "lwe/container/data_structure.hpp"
#include "lwe_params.hpp"
#include <stdexcept>

namespace LWE {

    template <typename T, typename Tp, class Params> class secret_key {
    public:
        Matrix<T, Params::pt_dim + Params::tau, Params::n> S_T;
        Matrix<Tp, Params::tau, Params::pt_dim> T_mat;

        secret_key() : S_T(), T_mat() {
            this->S_T.discrete_gaussian();
            this->T_mat.random_element();
        }
        secret_key(const secret_key &) = default;
        secret_key(secret_key &&) noexcept = default;
        secret_key &operator=(const secret_key &) = default;
        secret_key &operator=(secret_key &&) noexcept = default;
    };

    template <typename T, class Params> class public_parameter {
    public:
        Matrix<T, Params::n, Params::n> A;
        Matrix<T, Params::pt_dim + Params::tau, Params::n> D;

        public_parameter() : A(), D() {
            this->A.random_element();
            this->D.discrete_gaussian();
            this->D *= T(Params::p_int);
        }
        public_parameter(const public_parameter &) = default;
        public_parameter(public_parameter &&) noexcept = default;
        public_parameter &operator=(const public_parameter &) = default;
        public_parameter &operator=(public_parameter &&) noexcept = default;
    };

    template <typename T, typename Tp, class Params> class ciphertext {
    public:
        Vector<T, Params::n> a_vec;
        Vector<T, Params::pt_dim + Params::tau> c_vec;

        ciphertext() : a_vec(), c_vec() {}

        ciphertext(const ciphertext &other)
            : a_vec(other.a_vec), c_vec(other.c_vec) {}

        ciphertext &operator=(const ciphertext &other) {
            this->c_vec = other.c_vec;
            this->a_vec = other.a_vec;
            return *this;
        }

        ciphertext operator+(const ciphertext &other) const {
            ciphertext what(*this);
            what += other;
            return what;
        }

        ciphertext &operator+=(const ciphertext &other) {
            this->a_vec += other.a_vec;
            this->c_vec += other.c_vec;
            return *this;
        }

        ciphertext operator*(const Tp &val) const {
            ciphertext what(*this);
            what *= val;
            return what;
        }

        ciphertext &operator*=(const Tp &val) {
            T val_lifted;
            val.lift_to(val_lifted);
            this->a_vec *= val_lifted;
            this->c_vec *= val_lifted;
            return *this;
        }

        ciphertext &lift_ring_multiply(const Tp &val) {
            this->a_vec *= val;
            this->c_vec *= val;
            return *this;
        }

        ciphertext lift_ring_multiply(const Tp &val) const {
            ciphertext what(*this);
            return what.lift_ring_multiply(val);
        }

        ciphertext &rescale() {
            this->a_vec.rescale(Params::rescale_q, Params::p_int);
            this->c_vec.rescale(Params::rescale_q, Params::p_int);
            return *this;
        }
    };

    template <typename T, typename Tp, class Params>
    inline std::pair<secret_key<T, Tp, Params>, public_parameter<T, Params>>
    keygen() {
        secret_key<T, Tp, Params> sk;
        public_parameter<T, Params> pp;
        pp.D += sk.S_T * pp.A;
        return {sk, pp};
    }

    template <typename T, typename Tp, class Params>
    inline ciphertext<T, Tp, Params>
    encrypt(const secret_key<T, Tp, Params> &sk,
            const Vector<Tp, Params::pt_dim> &pt,
            bool ciphertext_noise = true) {
        Vector<Tp, Params::pt_dim + Params::tau> uv;
        for (uint32_t i = 0; i < Params::pt_dim; i++)
            uv[i] = pt[i];
        Vector<Tp, Params::tau> Tv = sk.T_mat * pt;
        for (uint32_t i = 0; i < Params::tau; i++)
            uv[i + Params::pt_dim] = Tv[i];

        ciphertext<T, Tp, Params> ct;
        ct.a_vec.random_element();

        for (uint32_t i = 0; i < Params::pt_dim + Params::tau; i++)
            uv[i].lift_to(ct.c_vec[i]);

        ct.c_vec += sk.S_T * ct.a_vec;

        Vector<T, Params::pt_dim + Params::tau> ev;
        if (ciphertext_noise) {
            ev.discrete_gaussian();
            ct.c_vec += ev * Params::p_int;
        }

        return ct;
    }

    template <typename T, typename Tp, class Params>
    inline Vector<Tp, Params::pt_dim>
    decrypt(const secret_key<T, Tp, Params> &sk,
            const ciphertext<T, Tp, Params> &ct,
            const uint128_t modulus_q = Params::q_int) {
        Vector<T, Params::pt_dim + Params::tau> zv =
            ct.c_vec - sk.S_T * ct.a_vec;

#ifdef DEBUG_LWE_ENCRYPTION_SCHEME
        for (uint32_t i = 0; i < Params::pt_dim + Params::tau; i++)
            std::cout << zv[i] << std::endl;
        std::cout << std::endl;
#endif

        Vector<Tp, Params::pt_dim + Params::tau> uv;
        for (uint32_t i = 0; i < Params::pt_dim + Params::tau; i++)
            uv[i].project_from(zv[i], modulus_q);

        Vector<Tp, Params::pt_dim> vv;
        for (uint32_t i = 0; i < Params::pt_dim; i++) {
            vv[i] = uv[i];
#ifdef DEBUG_LWE_ENCRYPTION_SCHEME
            std::cout << vv[i] << std::endl;
#endif
        }

        Vector<Tp, Params::tau> Tv;
        for (uint32_t i = 0; i < Params::tau; i++)
            Tv[i] = uv[Params::pt_dim + i];

        Vector<Tp, Params::tau> temp = sk.T_mat * vv;
        for (uint32_t i = 0; i < Params::tau; i++)
            if (temp[i] != Tv[i])
#ifdef DEBUG_LWE_ENCRYPTION_SCHEME
                std::cout << temp[i] << " " << Tv[i] << std::endl;
#else
                throw std::runtime_error("Decryption Corruption");
#endif

        return vv;
    }

    template <typename T, typename Tp, class Params>
    inline ciphertext<T, Tp, Params>
    add(const std::vector<Tp> &coefficients,
        const std::vector<ciphertext<T, Tp, Params>> &cts) {
        assert(cts.size() == coefficients.size());

        ciphertext<T, Tp, Params> res_ct;
        for (size_t i = 0; i < cts.size(); i++)
            res_ct += cts[i] * coefficients[i];

        return res_ct;
    }

    template <typename T, typename Tp, class Params>
    inline ciphertext<T, Tp, Params> &
    re_randomize(const public_parameter<T, Params> &pp,
                 ciphertext<T, Tp, Params> &ct) {
        Vector<T, Params::n> rv;
        rv.discrete_gaussian();

        ct.a_vec += pp.A * rv;
        ct.c_vec += pp.D * rv;

        Vector<T, Params::n> ea;
        ea.discrete_gaussian();
        ct.a_vec += ea * Params::p_int;

        Vector<T, Params::pt_dim + Params::tau> smudging_component;
        smudging_component.pm_bounded(Params::b_int);

        ct.c_vec += smudging_component * Params::p_int;

        return ct;
    }
}

#endif
