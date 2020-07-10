#ifndef __R1CS_LATTICE_SNARK__
#define __R1CS_LATTICE_SNARK__

#include "r1cs_lattice_snark_common.hpp"

namespace libsnark {

    /* DATA STRUCTURE DEFINITIONS */

    template <typename ppT, typename cpT, class Params>
    class r1cs_lattice_snark_crs {
    public:
        r1cs_constraint_system<libff::Fr<ppT>> constraint_system;
        std::vector<LWE::Vector<Rq_T<cpT>, Params::pt_dim + Params::tau>>
            enc_qs;
        LWE::public_parameter<Rq_T<cpT>, Params> public_parameter;
        LWERandomness::AES_KEY crs_aes_key{};

        r1cs_lattice_snark_crs() = default;
        r1cs_lattice_snark_crs &
        operator=(const r1cs_lattice_snark_crs &) = default;
        r1cs_lattice_snark_crs(const r1cs_lattice_snark_crs &) = default;
        r1cs_lattice_snark_crs(r1cs_lattice_snark_crs &&) noexcept = default;

        explicit r1cs_lattice_snark_crs(
            const r1cs_constraint_system<libff::Fr<ppT>> &cs,
            std::vector<LWE::Vector<Rq_T<cpT>, Params::pt_dim + Params::tau>>
                &&enc_q,
            LWE::public_parameter<Rq_T<cpT>, Params> &&pp,
            const LWERandomness::AES_KEY &_key)
            : constraint_system(cs), enc_qs(std::move(enc_q)),
              public_parameter(std::move(pp)), crs_aes_key{} {
            std::copy_n(_key.rd_key, 15, this->crs_aes_key.rd_key);
        }
    };

    template <typename ppT, typename cpT, class Params>
    class r1cs_lattice_snark_verification_key {
    public:
        LWE::secret_key<Rq_T<cpT>, libff::Fr<ppT>, Params> sk;
        std::vector<libff::Fr_vector<ppT>> A_prefix, B_prefix, C_prefix;
        libff::Fr_vector<ppT> Z_s;

        r1cs_lattice_snark_verification_key() = default;
        r1cs_lattice_snark_verification_key(
            LWE::secret_key<Rq_T<cpT>, libff::Fr<ppT>, Params> &&sk_,
            std::vector<libff::Fr_vector<ppT>> &&A_prefix_,
            std::vector<libff::Fr_vector<ppT>> &&B_prefix_,
            std::vector<libff::Fr_vector<ppT>> &&C_prefix_,
            libff::Fr_vector<ppT> &&Z_s_)
            : sk(std::move(sk_)), A_prefix(std::move(A_prefix_)),
              B_prefix(std::move(B_prefix_)), C_prefix(std::move(C_prefix_)),
              Z_s(std::move(Z_s_)) {}
    };

    template <typename ppT, typename Params>
    inline void
    gen_q_mat(const r1cs_constraint_system<libff::Fr<ppT>> &cs,
              r1cs_lattice_snark_query_matrix<ppT, Params::pt_dim> &q_mat,
              std::vector<libff::Fr_vector<ppT>> &A_qs,
              std::vector<libff::Fr_vector<ppT>> &B_qs,
              std::vector<libff::Fr_vector<ppT>> &C_qs,
              std::vector<libff::Fr_vector<ppT>> &A_prefix,
              std::vector<libff::Fr_vector<ppT>> &B_prefix,
              std::vector<libff::Fr_vector<ppT>> &C_prefix,
              std::vector<libff::Fr_vector<ppT>> &H_qs,
              libff::Fr_vector<ppT> &Z_s) {

        A_qs.resize(Params::query_num);
        B_qs.resize(Params::query_num);
        C_qs.resize(Params::query_num);
        H_qs.resize(Params ::query_num);

        Z_s.resize(Params::query_num);
        A_prefix.resize(Params::query_num);
        B_prefix.resize(Params::query_num);
        C_prefix.resize(Params::query_num);

        libff::enter_block("Generating QAP queries");
        size_t num_inputs = 0;

        auto t_s = reject_sampling_S(cs, Params::query_num);
        for (size_t i = 0; i < Params::query_num; i++) {
            qap_instance_evaluation<libff::Fr<ppT>> qap_inst =
                r1cs_to_qap_instance_map_with_evaluation(cs, t_s[i]);
            if (i == 0) {
                libff::print_indent();
                printf("* QAP number of variables: %zu\n",
                       qap_inst.num_variables());
                libff::print_indent();
                printf("* QAP pre degree: %zu\n", cs.constraints.size());
                libff::print_indent();
                printf("* QAP degree: %zu\n", qap_inst.degree());
                libff::print_indent();
                printf("* QAP number of input variables: %zu\n",
                       qap_inst.num_inputs());

                num_inputs = qap_inst.num_inputs();
            }
            A_qs[i] = std::move(qap_inst.At);
            B_qs[i] = std::move(qap_inst.Bt);
            C_qs[i] = std::move(qap_inst.Ct);
            H_qs[i] = std::move(qap_inst.Ht);
            Z_s[i] = qap_inst.Zt;
            A_prefix[i].reserve(num_inputs + 1);
            std::copy_n(std::begin(A_qs[i]), num_inputs + 1,
                        std::begin(A_prefix[i]));
            B_prefix[i].reserve(num_inputs + 1);
            std::copy_n(std::begin(B_qs[i]), num_inputs + 1,
                        std::begin(B_prefix[i]));
            C_prefix[i].reserve(num_inputs + 1);
            std::copy_n(std::begin(C_qs[i]), num_inputs + 1,
                        std::begin(C_prefix[i]));
        }

        const uint64_t ABC_rows = A_qs.begin()->size() - num_inputs - 1;
        const uint64_t H_rows = H_qs.begin()->size();
        const uint64_t rows = ABC_rows + 3 + H_rows;
        q_mat.resize(rows);
        for (uint64_t i = 0; i < ABC_rows; i++) {
            for (uint64_t j = 0; j < Params::query_num; j++) {
                q_mat[i][j * LWE::query_size] = A_qs[j][i + 1 + num_inputs];
                q_mat[i][j * LWE::query_size + 1] = B_qs[j][i + 1 + num_inputs];
                q_mat[i][j * LWE::query_size + 2] = C_qs[j][i + 1 + num_inputs];
            }
        }
        for (uint64_t i = 0; i < 3; i++)
            for (uint64_t j = 0; j < Params::query_num; j++)
                q_mat[ABC_rows + i][i + LWE::query_size * j] = Z_s[j];
        for (uint64_t i = 0; i < Params::query_num; i++)
            for (uint64_t j = 0; j < H_rows; j++)
                q_mat[ABC_rows + 3 + j][i * LWE::query_size + 3] = H_qs[i][j];
        libff::leave_block("Generating QAP queries");
    }

    template <typename ppT, typename cpT, class Params>
    void r1cs_lattice_snark_generator(
        const r1cs_constraint_system<libff::Fr<ppT>> &cs,
        r1cs_lattice_snark_crs<ppT, cpT, Params> &crs,
        r1cs_lattice_snark_verification_key<ppT, cpT, Params> &vk) {

        libff::enter_block("Generating LWE secret key");
        auto sk_pp = LWE::keygen<Rq_T<cpT>, libff::Fr<ppT>, Params>();
        libff::leave_block("Generating LWE secret key");

        std::vector<libff::Fr_vector<ppT>> A_qs, B_qs, C_qs, H_qs;
        libff::Fr_vector<ppT> Zs;
        std::vector<libff::Fr_vector<ppT>> A_prefix, B_prefix, C_prefix;
        r1cs_lattice_snark_query_matrix<ppT, Params::pt_dim> q_mat;
        gen_q_mat<ppT, Params>(cs, q_mat, A_qs, B_qs, C_qs, A_prefix, B_prefix,
                               C_prefix, H_qs, Zs);

        libff::enter_block("Generating CRS and VK");
        LWERandomness::AES_KEY _crs_aes_key;
        genAES_key(&_crs_aes_key);

        std::vector<LWE::Vector<Rq_T<cpT>, Params::pt_dim + Params::tau>> dummy;
        crs = r1cs_lattice_snark_crs<ppT, cpT, Params>(
            cs, std::move(dummy), std::move(sk_pp.second), _crs_aes_key);
        vk = r1cs_lattice_snark_verification_key<ppT, cpT, Params>(
            std::move(sk_pp.first), std::move(A_prefix), std::move(B_prefix),
            std::move(C_prefix), std::move(Zs));
        encrypt_query_matrix<ppT, cpT, Params>(vk.sk, q_mat, crs.crs_aes_key,
                                               crs.enc_qs);

        libff::leave_block("Generating CRS and VK");
    }

    template <typename ppT>
    inline void prepare_pi_proof(const qap_witness<libff::Fr<ppT>> &qap_wit,
                                 libff::Fr_vector<ppT> &pi) {
        libff::enter_block("Prepare pi proof");
        size_t num_inputs = qap_wit.num_inputs();
        size_t num_ABC_coeffs =
            qap_wit.coefficients_for_ABCs.size() - num_inputs;
        size_t proof_dim =
            num_ABC_coeffs + 3 + qap_wit.coefficients_for_H.size();
        pi.resize(proof_dim);

        for (size_t i = 0; i < num_ABC_coeffs; i++)
            pi[i] = qap_wit.coefficients_for_ABCs[i + num_inputs];
        pi[num_ABC_coeffs] = qap_wit.d1;
        pi[num_ABC_coeffs + 1] = qap_wit.d2;
        pi[num_ABC_coeffs + 2] = qap_wit.d3;
        std::copy(std::begin(qap_wit.coefficients_for_H),
                  std::end(qap_wit.coefficients_for_H),
                  std::begin(pi) + num_ABC_coeffs + 3);
        libff::leave_block("Prepare pi proof");
    }

    template <typename ppT, typename cpT, class Params>
    r1cs_lattice_snark_proof<ppT, cpT, Params> r1cs_lattice_snark_prove(
        const r1cs_lattice_snark_crs<ppT, cpT, Params> &crs,
        const r1cs_primary_input<libff::Fr<ppT>> &primary_input,
        const r1cs_auxiliary_input<libff::Fr<ppT>> &auxiliary_input) {
        libff::enter_block("Call to r1cs lattice snark prover");

        libff::enter_block("Compute H polynomial");

        LWERandomness::AES_KEY _aes_key;
        genAES_key(&_aes_key);

        auto *temp_prg = new LWERandomness::PseudoRandomGenerator(_aes_key);
        auto *temp_dg = new LWERandomness::DiscreteGaussian(
            Params::width, LWE::expand, *temp_prg);
        auto *original_prg = ppT::prg;
        auto *original_dg = ppT::dg;
        public_params_init<ppT, cpT>(temp_prg, temp_dg);

        const libff::Fr<ppT> d1 = libff::Fr<ppT>::random_element(),
                             d2 = libff::Fr<ppT>::random_element(),
                             d3 = libff::Fr<ppT>::random_element();
        public_params_init<ppT, cpT>(original_prg, original_dg);
        delete temp_dg;
        delete temp_prg;

        const qap_witness<libff::Fr<ppT>> qap_wit = r1cs_to_qap_witness_map(
            crs.constraint_system, primary_input, auxiliary_input, d1, d2, d3);
        libff::leave_block("Compute H polynomial");

        libff::Fr_vector<ppT> pi;
        prepare_pi_proof<ppT>(qap_wit, pi);
        assert(pi.size() == crs.enc_qs->size());

        temp_prg = new LWERandomness::PseudoRandomGenerator(crs.crs_aes_key);
        temp_dg = new LWERandomness::DiscreteGaussian(Params::width,
                                                      LWE::expand, *temp_prg);
        original_prg = ppT::prg;
        original_dg = ppT::dg;
        public_params_init<ppT, cpT>(temp_prg, temp_dg);

        libff::enter_block("Generating response");

        LWE::ciphertext<Rq_T<cpT>, libff::Fr<ppT>, Params> added;

        int index = crs.enc_qs.size();

        libff::enter_block("Generating ct linear combination");
        for (int i = 0; i < index; i++)
            added.c_vec += crs.enc_qs[i] * pi[i];
        libff::leave_block("Generating ct linear combination");

        LWE::Vector<Rq_T<cpT>, Params::n> a_vec;

        std::chrono::duration<double, std::micro> dur_rand{};
        std::chrono::duration<double, std::milli> dur_lin{};
        libff::enter_block("Generating a_vec linear combination");

        for (int i = 0; i < index; i++) {
            auto rnd_srt = std::chrono::high_resolution_clock::now();
            a_vec.random_element();
            auto rnd_end = std::chrono::high_resolution_clock::now();
            auto linear_srt = std::chrono::high_resolution_clock::now();
            added.a_vec += (a_vec *= pi[i]);
            auto linear_end = std::chrono::high_resolution_clock::now();
            dur_rand += (rnd_end - rnd_srt);
            dur_lin += (linear_end - linear_srt);
        }

        using milli_s = std::chrono::milliseconds;
        double rand_t = std::chrono::duration_cast<milli_s>(dur_rand).count(),
               linear_t = std::chrono::duration_cast<milli_s>(dur_lin).count();
        libff::leave_block("Generating a_vec linear combination");

        public_params_init<ppT, cpT>(original_prg, original_dg);
        delete temp_dg;
        delete temp_prg;

#ifndef NOT_PROVABLE_ZK
        LWE::re_randomize(crs.public_parameter, added);
#endif
        added.rescale();
        libff::leave_block("Generating response");

        libff::leave_block("Call to r1cs lattice snark prover");
        std::cout << "\n  * Randomize a_vec " +
                         std::to_string((double) rand_t / 1000.0) + "s\n"
                  << "  * Linear comb on a_vec " +
                         std::to_string((double) linear_t / 1000.0) + "s\n"
                  << "  * Linear comb size " + std::to_string(index)
                  << std::endl;
        return r1cs_lattice_snark_proof<ppT, cpT, Params>(std::move(added));
    }

    template <typename ppT, typename cpT, class Params>
    bool r1cs_lattice_snark_verify(
        const r1cs_lattice_snark_verification_key<ppT, cpT, Params> &vk,
        const r1cs_primary_input<libff::Fr<ppT>> &primary_input,
        const r1cs_lattice_snark_proof<ppT, cpT, Params> &proof) {
        bool res = true;
        libff::enter_block("Call to r1cs lattice snark verifier");

        libff::enter_block("Decrypting proof");
        auto decrypted = LWE::decrypt(vk.sk, proof.response, Params::rescale_q);
        libff::Fr_vector<ppT> Ap(Params::query_num), Bp(Params::query_num),
            Cp(Params::query_num), Hp(Params::query_num);

        for (uint i = 0; i < Params::query_num; i++) {
            Ap[i] = decrypted[i * LWE::query_size];
            Bp[i] = decrypted[i * LWE::query_size + 1];
            Cp[i] = decrypted[i * LWE::query_size + 2];
            Hp[i] = decrypted[i * LWE::query_size + 3];

            Ap[i] += vk.A_prefix[i][0];
            Bp[i] += vk.B_prefix[i][0];
            Cp[i] += vk.C_prefix[i][0];

            for (uint64_t j = 0; j < primary_input.size(); j++) {
                Ap[i] += primary_input[j] * vk.A_prefix[i][j + 1];
                Bp[i] += primary_input[j] * vk.B_prefix[i][j + 1];
                Cp[i] += primary_input[j] * vk.C_prefix[i][j + 1];
            }
        }
        libff::leave_block("Decrypting proof");

        libff::enter_block("Checking QAP divisibility");
        for (uint i = 0; i < Params::query_num; i++) {
            if (Ap[i] * Bp[i] != Hp[i] * vk.Z_s[i] + Cp[i]) {
                if (!libff::inhibit_profiling_info) {
                    libff::print_indent();
                    printf("QAP divisibility check failed.\n");
                }
                res = false;
            }
        }
        libff::leave_block("Checking QAP divisibility");

        libff::leave_block("Call to r1cs lattice snark verifier");
        return res;
    }

    template <typename ppT, typename cpT, class Params>
    bool run_r1cs_lattice_snark(const r1cs_example<libff::Fr<ppT>> &example) {
        libff::enter_block("Call to R1CS lattice SNARK");

        libff::print_header("R1CS lattice SNARK Generator");
        r1cs_lattice_snark_crs<ppT, cpT, Params> crs;
        r1cs_lattice_snark_verification_key<ppT, cpT, Params> vk;
        r1cs_lattice_snark_generator<ppT, cpT, Params>(
            example.constraint_system, crs, vk);
        printf("\n");
        libff::print_indent();
        libff::print_mem("after generator");

        libff::print_header("R1CS lattice SNARK Prover");
        r1cs_lattice_snark_proof<ppT, cpT, Params> proof =
            r1cs_lattice_snark_prove<ppT, cpT, Params>(
                crs, example.primary_input, example.auxiliary_input);
        printf("\n");
        libff::print_indent();
        libff::print_mem("after prover");

        libff::print_header("R1CS lattice SNARK Verifier");
        const bool ans = r1cs_lattice_snark_verify<ppT, cpT>(
            vk, example.primary_input, proof);
        printf("\n");
        libff::print_indent();
        libff::print_mem("after verifier");
        printf("* The verification result is: %s\n", (ans ? "PASS" : "FAIL"));

        libff::leave_block("Call to R1CS lattice SNARK");
        return ans;
    }

    template <typename ppT, typename cpT, class Params>
    void test_r1cs_lattice_snark(size_t num_constraints, size_t input_size) {
        libff::print_header("(enter) Test R1CS lattice SNARK");
        r1cs_example<libff::Fr<ppT>> example =
            generate_r1cs_example_with_field_input<libff::Fr<ppT>>(
                num_constraints, input_size);
        const bool res = run_r1cs_lattice_snark<ppT, cpT, Params>(example);
        if (!res)
            libff::print_header("TEST FAILED");

        libff::print_header("(leave) Test R1CS lattice SNARK");
    }
}

#endif
