add_library(
        lwe_interface
        STATIC

        lwe.hpp
        lwe_params.hpp
        snark/r1cs_lattice_snark_common.hpp
        snark/r1cs_lattice_snark.hpp

        tests/common.hpp

        container/data_structure.hpp
        container/ring_base.hpp
        container/extension.hpp
        container/field_base.hpp
        container/utils.hpp
        randomness/prg.hpp
        randomness/aes.hpp
)

target_link_libraries(lwe_interface snark)

target_include_directories(
        lwe_interface

        PUBLIC
        ${LIBSNARK_DIR}
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

add_executable(r1cs_lattice_snark_prime19_test tests/SNARK/test_r1cs_lattice_snark_p19.cpp)

target_link_libraries(r1cs_lattice_snark_prime19_test lwe_interface)

add_executable(r1cs_lattice_snark_lpcp_tr_prime19_test tests/SNARK/test_r1cs_lattice_snark_p19_lpcp_tr.cpp)

target_link_libraries(r1cs_lattice_snark_lpcp_tr_prime19_test lwe_interface)

add_executable(r1cs_lattice_snark_prime13_test tests/SNARK/test_r1cs_lattice_snark_p13.cpp)

target_link_libraries(r1cs_lattice_snark_prime13_test lwe_interface)

add_executable(r1cs_lattice_snark_lpcp_tr_prime13_test tests/SNARK/test_r1cs_lattice_snark_p13_lpcp_tr.cpp)

target_link_libraries(r1cs_lattice_snark_lpcp_tr_prime13_test lwe_interface)

add_executable(r1cs_lattice_snark_prime23_test tests/SNARK/test_r1cs_lattice_snark_p23.cpp)

target_link_libraries(r1cs_lattice_snark_prime23_test lwe_interface)

add_executable(r1cs_lattice_snark_prime28_test tests/SNARK/test_r1cs_lattice_snark_p28.cpp)

target_link_libraries(r1cs_lattice_snark_prime28_test lwe_interface)
