#ifndef __LWE_TEST_COMMON__
#define __LWE_TEST_COMMON__

#define CURVE_BN128

#include "lwe/container/extension.hpp"
#include "lwe/container/field_base.hpp"
#include "lwe/container/ring_base.hpp"
#include "lwe/lwe_params.hpp"
#include "lwe/tests/circ_lattice_params.hpp"

#include <libff/algebra/fields/fp.hpp>
#include <libff/common/profiling.hpp>
#include <libsnark/common/default_types/r1cs_gg_ppzksnark_pp.hpp>
#include <libsnark/reductions/r1cs_to_qap/r1cs_to_qap.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/examples/run_r1cs_gg_ppzksnark.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>

class Fp2_b13_pp {
public:
    using Fp = libsnark::Field<uint32_t, LWE::B13Fp2ParamsBase::p_int>;
    using Fp_type = libsnark::Extension<Fp>;

    static LWERandomness::PseudoRandomGenerator *prg;
    static LWERandomness::DiscreteGaussian *dg;

    static void init_public_params() {
        Fp::prg = Fp2_b13_pp::prg;
        Fp::dg = Fp2_b13_pp::dg;
        Fp::multiplicative_generator = Fp(3);
        Fp::s = 1;
        Fp::root_of_unity = Fp(8190);

        Fp_type::non_residue = Fp(3);
        Fp_type::multiplicative_generator = Fp_type(3, 1);
        Fp_type::s = 14;
        Fp_type::root_of_unity = Fp_type(8127, 64);
    }
};

LWERandomness::PseudoRandomGenerator *Fp2_b13_pp::prg;
LWERandomness::DiscreteGaussian *Fp2_b13_pp::dg;

class Fp2_b19_pp {
public:
    using Fp = libsnark::Field<uint64_t, LWE::B19Fp2ParamsBase::p_int>;
    using Fp_type = libsnark::Extension<Fp>;

    static LWERandomness::PseudoRandomGenerator *prg;
    static LWERandomness::DiscreteGaussian *dg;

    static void init_public_params() {
        Fp::prg = Fp2_b19_pp::prg;
        Fp::dg = Fp2_b19_pp::dg;
        Fp::multiplicative_generator = Fp(3);
        Fp::s = 1;
        Fp::root_of_unity = Fp(524286);

        Fp_type::non_residue = Fp(3);
        Fp_type::multiplicative_generator = Fp_type(3, 1);
        Fp_type::s = 20;
        Fp_type::root_of_unity = Fp_type(512, 523775);
    }
};

LWERandomness::PseudoRandomGenerator *Fp2_b19_pp::prg;
LWERandomness::DiscreteGaussian *Fp2_b19_pp::dg;

template <LWE::uint128_t q_int = LWE::B19C20::q_int> class Ring2_common_pp {
public:
    using Rq = libsnark::Ring<LWE::uint128_t, q_int>;
    using Rq_type = libsnark::Extension<Rq>;

    static LWERandomness::PseudoRandomGenerator *prg;
    static LWERandomness::DiscreteGaussian *dg;

    static void init_public_params() {
        Rq::prg = prg;
        Rq::dg = dg;

        Rq_type::non_residue = Rq(3);
    }
};

template <LWE::uint128_t q_int>
LWERandomness::PseudoRandomGenerator *Ring2_common_pp<q_int>::prg;
template <LWE::uint128_t q_int>
LWERandomness::DiscreteGaussian *Ring2_common_pp<q_int>::dg;

class Fp_b28_pp {
public:
    using Fp_type = libsnark::Field<uint64_t, LWE::B28FpParamsBase::p_int>;

    static LWERandomness::PseudoRandomGenerator *prg;
    static LWERandomness::DiscreteGaussian *dg;

    static void init_public_params() {
        Fp_type::prg = Fp_b28_pp::prg;
        Fp_type::dg = Fp_b28_pp::dg;
        Fp_type::s = 25;
        Fp_type::multiplicative_generator = Fp_type(3);
        Fp_type::root_of_unity = Fp_type(243);
    }
};

LWERandomness::PseudoRandomGenerator *Fp_b28_pp::prg;
LWERandomness::DiscreteGaussian *Fp_b28_pp::dg;

class Fp_b23_pp {
public:
    using Fp_type = libsnark::Field<uint64_t, LWE::B23FpParamsBase::p_int>;

    static LWERandomness::PseudoRandomGenerator *prg;
    static LWERandomness::DiscreteGaussian *dg;

    static void init_public_params() {
        Fp_type::prg = Fp_b23_pp::prg;
        Fp_type::dg = Fp_b23_pp::dg;
        Fp_type::s = 20;
        Fp_type::multiplicative_generator = Fp_type(3);
        Fp_type::root_of_unity = Fp_type(2187);
    }
};

LWERandomness::PseudoRandomGenerator *Fp_b23_pp::prg;
LWERandomness::DiscreteGaussian *Fp_b23_pp::dg;

class Fp_b19_pp {
public:
    using Fp_type = libsnark::Field<uint64_t, LWE::B19FpParamsBase::p_int>;

    static LWERandomness::PseudoRandomGenerator *prg;
    static LWERandomness::DiscreteGaussian *dg;

    static void init_public_params() {
        Fp_type::prg = Fp_b19_pp::prg;
        Fp_type::dg = Fp_b19_pp::dg;
        Fp_type::s = 1;
        Fp_type::multiplicative_generator = Fp_type(3);
        Fp_type::root_of_unity = Fp_type(524286);
    }
};

LWERandomness::PseudoRandomGenerator *Fp_b19_pp::prg;
LWERandomness::DiscreteGaussian *Fp_b19_pp::dg;

class Fp_b13_pp {
public:
    using Fp_type = libsnark::Field<uint32_t, LWE::B13FpParamsBase::p_int>;

    static LWERandomness::PseudoRandomGenerator *prg;
    static LWERandomness::DiscreteGaussian *dg;

    static void init_public_params() {
        Fp_type::prg = Fp_b13_pp::prg;
        Fp_type::dg = Fp_b13_pp::dg;
        Fp_type::multiplicative_generator = Fp_type(3);
        Fp_type::s = 1;
        Fp_type::root_of_unity = Fp_type(8190);
    }
};

LWERandomness::PseudoRandomGenerator *Fp_b13_pp::prg;
LWERandomness::DiscreteGaussian *Fp_b13_pp::dg;

template <LWE::uint128_t q_int = LWE::B19C20::q_int> class Ring_common_pp {
public:
    using Rq_type = libsnark::Ring<LWE::uint128_t, q_int>;
    using Rq = Rq_type;

    static LWERandomness::PseudoRandomGenerator *prg;
    static LWERandomness::DiscreteGaussian *dg;

    static void init_public_params() {
        Rq_type::prg = prg;
        Rq_type::dg = dg;
    }
};

template <LWE::uint128_t q_int>
LWERandomness::PseudoRandomGenerator *Ring_common_pp<q_int>::prg;
template <LWE::uint128_t q_int>
LWERandomness::DiscreteGaussian *Ring_common_pp<q_int>::dg;

namespace LWE {
    constexpr const double width_default = 64.0;
}

#endif
