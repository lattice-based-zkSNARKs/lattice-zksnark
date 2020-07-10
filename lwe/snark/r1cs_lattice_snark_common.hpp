#ifndef __R1CS_LATTICE_SNARK_COMMON__
#define __R1CS_LATTICE_SNARK_COMMON__

#include "lwe/container/extension.hpp"
#include "lwe/lwe.hpp"
#include "lwe/lwe_params.hpp"
#include "lwe/randomness/aes.hpp"
#include "lwe/randomness/prg.hpp"

#include <chrono>
#include <libff/algebra/curves/public_params.hpp>
#include <libfqfft/evaluation_domain/get_evaluation_domain.hpp>
#include <libsnark/reductions/r1cs_to_qap/r1cs_to_qap.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/r1cs.hpp>

namespace libsnark {

    /* TYPE ALIAS DEFINITONS */

    template <typename cpT> using Rq_T = typename cpT::Rq_type;

    template <typename ppT, uint32_t pt_dim>
    using r1cs_lattice_snark_query_matrix =
        std::vector<LWE::Vector<libff::Fr<ppT>, pt_dim>>;

    template <typename ppT, typename cpT, class Params>
    class r1cs_lattice_snark_proof {
    public:
        LWE::ciphertext<Rq_T<cpT>, libff::Fr<ppT>, Params> response;
        r1cs_lattice_snark_proof() = default;
        explicit r1cs_lattice_snark_proof(
            LWE::ciphertext<Rq_T<cpT>, libff::Fr<ppT>, Params> &&response)
            : response(std::move(response)) {}
    };

    template <typename ppT>
    std::vector<ppT> reject_sampling_S(const r1cs_constraint_system<ppT> &cs,
                                       int sample_num) {
        const auto domain = libfqfft::get_evaluation_domain<ppT>(
            cs.num_constraints() + cs.num_inputs() + 1);
        std::vector<ppT> res(sample_num);
        for (int i = 0; i < sample_num; i++) {
            ppT _res = ppT::random_element();
            while (domain->compute_vanishing_polynomial(_res) == ppT::zero())
                _res = ppT::random_element();
            res[i] = _res;
        }
        return res;
    }

    template <typename T, typename... RT>
    inline void public_params_init(LWERandomness::PseudoRandomGenerator *prg,
                                   LWERandomness::DiscreteGaussian *dg) {
        T::prg = prg;
        T::dg = dg;
        T::init_public_params();
        if constexpr (sizeof...(RT) > 0)
            public_params_init<RT...>(prg, dg);
    }

    inline void genAES_key(LWERandomness::AES_KEY *_key) {
        static std::ifstream urandom("/dev/urandom", std::ios::binary);
        LWERandomness::byte buffer[LWERandomness::AES_KEY_BYTES];
        urandom.read(reinterpret_cast<char *>(buffer),
                     LWERandomness::AES_KEY_BYTES);
        urandom.close();
        LWERandomness::AES_128_Key_Expansion(buffer, _key);
    }

    template <typename ppT, typename cpT, class Params>
    static inline void encrypt_query_matrix(
        const LWE::secret_key<Rq_T<cpT>, libff::Fr<ppT>, Params> &sk,
        const r1cs_lattice_snark_query_matrix<ppT, Params::pt_dim> &q_mat,
        const LWERandomness::AES_KEY &crs_aes_key,
        std::vector<LWE::Vector<Rq_T<cpT>, Params::pt_dim + Params::tau>>
            &enc_qs) {
        enc_qs.resize(q_mat.size());
        auto *temp_prg = new LWERandomness::PseudoRandomGenerator(crs_aes_key);
        auto *temp_dg = new LWERandomness::DiscreteGaussian(
            Params::width, LWE::expand, *temp_prg);
        auto *original_prg = ppT::prg;
        auto *original_dg = ppT::dg;
        public_params_init<ppT, cpT>(temp_prg, temp_dg);

        int counter = 0;
        for (const auto &row : q_mat) {
            auto encrypted_query =
                LWE::encrypt<Rq_T<cpT>, libff::Fr<ppT>, Params>(sk, row, false);
            enc_qs[counter++] = std::move(encrypted_query.c_vec);
        }

        public_params_init<ppT, cpT>(original_prg, original_dg);
        for (auto &row : enc_qs) {
            LWE::Vector<Rq_T<cpT>, Params::pt_dim + Params::tau> ev;
            ev.discrete_gaussian();
            row += ev * Params::p_int;
        }

        delete temp_dg;
        delete temp_prg;
    }

    template <typename pdpT, typename ppT, uint32_t pt_dim>
    inline void expand_queries(
        const r1cs_lattice_snark_query_matrix<pdpT, pt_dim> &orig_query,
        r1cs_lattice_snark_query_matrix<ppT, pt_dim * 2> &expand_query) {
        size_t orig_query_len = orig_query.size();
        expand_query.resize(orig_query_len * 2);
        for (size_t i = 0; i < orig_query_len; i++) {
            for (size_t j = 0; j < pt_dim; j++) {
                expand_query[2 * i][2 * j] = orig_query[i][j].c0;
                expand_query[2 * i][2 * j + 1] = orig_query[i][j].c1;
                expand_query[2 * i + 1][2 * j] =
                    orig_query[i][j].c1 * libff::Fr<pdpT>::non_residue;
                expand_query[2 * i + 1][2 * j + 1] = orig_query[i][j].c0;
            }
        }
    }

    template <typename ppT>
    inline void fp_shrink(const std::vector<ppT> &fp_respond,
                          std::vector<libsnark::Extension<ppT>> &shrink_res) {
        shrink_res.resize(fp_respond.size() / 2);
        for (size_t i = 0; i < shrink_res.size(); i++) {
            shrink_res[i].c0 = fp_respond[i * 2];
            shrink_res[i].c1 = fp_respond[i * 2 + 1];
        }
    }
}

#endif
