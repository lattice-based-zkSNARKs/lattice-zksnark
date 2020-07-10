# Lattice-Based zkSNARKs over libsnark

**WARNING: This code provides a proof-of-concept implementation of the lattice-based zero-knowledge succinct non-interactive argument of knowledge (zkSNARK). This is a research prototype and is not intended to be used directly in critical or production-level systems.**

This is the implementation of the construction described in the following CCS 2021 paper: [**Shorter and Faster Post-Quantum Designated-Verifier zkSNARKs from Lattices**](https://eprint.iacr.org/2021/977).

The implementation is built on the top of the [libsnark](https://github.com/ISW21-Implementation/libsnark) library. All of our code is released under the MIT License (see [LICENSE](./LICENSE)).

This implementation is currently maintained by Hang Su (hs2nu@virginia.edu) and David Wu (dwu4@cs.utexas.edu).

---

## Overview

We refer to [libsnark](https://github.com/lattice-based-zkSNARKs/libsnark) for an overview of zero-knowledge succinct non-interactive arguments of knowledge (zkSNARKs). This repository provides an implementation of the lattice-based zkSNARK described in [ISW21] for the language of rank-1 constraint satisfiability (R1CS).

The parameters in this prototype implementation support R1CS instances with roughly `2^20` constraints. We set the computational security parameter to 128 and the statistical security parameter (for zero knowledge) to 40. The parameters are based on the best quantum attacks on lattice-based cryptography (see [LWE Estimator](https://lwe-estimator.readthedocs.io/en/latest/#)).

---

## Build Instructions

### Dependencies

The lattice-based zkSNARKs library relies on the following:

- C++ build environment supporting intrinsic `__int128_t` and `__uint128_t`
- CMake build infrastructure
- Git
- libsnark (fetched via Git submodules)
- Processor supporting AES-NI instruction set

### Building

Fetch the libsnark submodule by

```shell
git submodule update --init --recursive
```

Create the Makefile by

```shell
mkdir build && cd build && cmake ..
```

Compile the library and proceed the tests, run `make` command in the `build` directory.

---

## Running Example

Under `./build`, run

```shell
./lwe/bin/r1cs_lattice_snark_prime19_test 10000 100
```

The command tests the zkSNARK (over the quadratic extension field with characteristic `p = 2^19 - 1`) by

- sampling a CRS and a secret verification key for verifying R1CS instances with 10000 constraints and statements of length 100
- invoking the prover on a statement and the verifier on the resulting proof

The R1CS instance is drawn from the example in libsnark.

---

## Construction Settings

This library provides implementations of several instantiations considered in [ISW21]. At a high-level,
the [ISW21] construction follows the [BCIOP13, GGPR13, BISW17] approach of compiling a linear probabilistically
checkable proof (linear PCP) into a zkSNARK using linear-only vector encryption. The underlying linear PCP
is based on quadratic arithmetic programs [GGPR13].

### Linear PCPs and Linear-Only Vector Encryption over Quadratic Extensions

The main instantiations consider a linear PCP and a linear-only vector encryption over quadratic extension
fields. This library supports quadratic extensions with characteristic 2^13 - 1 and 2^19 - 1. These
constructions achieve the smallest proofs and best performance. They rely on module lattices (of rank 2).

| Field Characteristic | Executable (inside `./build/`) |
| -------------------- | ------------------------------ |
| [`2^13 - 1`](lwe/tests/SNARK/test_r1cs_lattice_snark_p13.cpp) | `./lwe/bin/r1cs_lattice_snark_prime13_test`         |
| [`2^19 - 1`](lwe/tests/SNARK/test_r1cs_lattice_snark_p19.cpp) | `./lwe/bin/r1cs_lattice_snark_prime19_test`         |

### Linear PCP over Quadratic Extension Field, Linear-Only Vector Encryption over Base Field

An alternative instantiation considers a linear PCP over a quadratic extension field and
a linear-only vector encryption over the base field. These constructions first apply a generic transformation
to the linear PCP over the extension field to obtain a (longer) linear PCP over the base field. The advantage
of these constructions is that they only require integer lattices rather than module lattices:

| Field Characteristic | Executable (inside `./build/`) |
| -------------------- | ------------------------------ |
| [`2^13 - 1`](lwe/tests/SNARK/test_r1cs_lattice_snark_p13_lpcp_tr.cpp) | `./lwe/bin/r1cs_lattice_snark_lpcp_tr_prime13_test`   |
| [`2^19 - 1`](lwe/tests/SNARK/test_r1cs_lattice_snark_p19_lpcp_tr.cpp) | `./lwe/bin/r1cs_lattice_snark_lpcp_tr_prime19_test` |

### Linear PCP and Linear-Only Vector Encryption over Base Field

The final instantiations we support directly compile linear PCPs over a (larger) base field with a linear-only
vector encryption over the same base field. These instantiations also work over integer lattices, but
have longer proofs and in general, worse efficiency compared to the constructions over extension fields.

| Field Characteristic | Executable (inside `./build/`) |
| -------------------- | ------------------------------ |
| [`7 x 2^20 + 1`](lwe/tests/SNARK/test_r1cs_lattice_snark_p23.cpp) | `./lwe/bin/r1cs_lattice_snark_prime23_test`         |
| [`5 x 2^25 + 1`](lwe/tests/SNARK/test_r1cs_lattice_snark_p28.cpp) | `./lwe/bin/r1cs_lattice_snark_prime28_test`         |

All of the testing executables follow the input format: `executable_file [constraint_size] [input_size]`.

The parameters for the above schemes are chosen to support R1CS systems with roughly 2^20 constraints.

---

## Citation

```Bibtex
@inproceedings{ISW21,
  author    = {Yuval Ishai and Hang Su and David J. Wu},
  title     = {Shorter and Faster Post-Quantum Designated-Verifier {zkSNARKs} from Lattices},
  booktitle = {{ACM} {CCS}},
  year      = {2021}
}
```

---

## References

\[BCIOP13\] [Succinct Non-interactive Arguments via Linear Interactive Proofs](https://eprint.iacr.org/2012/718.pdf). Nir Bitansky, Alessandro Chiesa, Yuval Ishai, Rafail Ostrovsky, and Omer Paneth. TCC, 2013.

\[BISW17\] [Lattice-Based SNARGs and Their Application to More Efficient Obfuscation](https://eprint.iacr.org/2017/240). Dan Boneh, Yuval Ishai, Amit Sahai, and David J. Wu. EUROCRYPT, 2017.

\[GGPR13\] [Quadratic Span Programs and Succinct NIZKs without PCPs](https://eprint.iacr.org/2012/215.pdf). Rosario Gennaro, Craig Gentry, Bryan Parno, and Mariana Raykova. EUROCRYPT, 2013.

\[ISW21\] [Shorter and Faster Post-Quantum Designated-Verifier zkSNARKs from Lattices](https://eprint.iacr.org/2021/977.pdf). Yuval Ishai, Hang Su, and David J. Wu. ACM CCS, 2021.

\[PVW08\] [A Framework for Efficient and Composable Oblivious Transfer](https://eprint.iacr.org/2007/348.pdf). Chris Peikert, Vinod Vaikuntanathan, and Brent Waters. CRYPTO, 2008.
