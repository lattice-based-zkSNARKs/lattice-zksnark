#include "lwe/snark/r1cs_lattice_snark.hpp"
#include "lwe/tests/circ_lattice_params.hpp"
#include "lwe/tests/common.hpp"

using namespace libsnark;
using namespace LWE;
using Param = LWE::B13C20;

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cout << "usage: ./lwe/bin/r1cs_lattice_snark constraint_n input_n"
                  << std::endl;
        return -1;
    }

    using ring_pp = Ring2_common_pp<Param::q_int>;

    auto prg = new LWERandomness::PseudoRandomGenerator();
    auto dg =
        new LWERandomness::DiscreteGaussian(Param::width, LWE::expand, *prg);

    public_params_init<Fp2_b13_pp, ring_pp>(prg, dg);

    int num_constraints = std::stoi(std::string(argv[1])),
        input_size = std::stoi(std::string(argv[2]));

    std::cout << "q_log = " << Param::q_log << "\n"
              << "q_rescale = " << Param::rescale_q << "\n"
              << "b = " << Param::b_int << "\n"
              << "pt_dim = " << Param::pt_dim << "\n"
              << "n = " << Param::n << "\n"
              << "s = " << Param::width << std::endl;

    libff::start_profiling();
    test_r1cs_lattice_snark<Fp2_b13_pp, ring_pp, Param>(num_constraints,
                                                        input_size);

    delete prg;
    delete dg;

    return 0;
}
