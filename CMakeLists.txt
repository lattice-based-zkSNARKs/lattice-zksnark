cmake_minimum_required(VERSION 3.11)
project (lwe_interface)

if("${LIBSNARK_DIR}" STREQUAL "")
  set(
    LIBSNARK_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/depends/libsnark"
    CACHE
    STRING
    "Optionally specify the libsnark installation directory relative to the source directory (default: inside libsnark folder)"
	)
else()
  set(LIBSNARK_DIR "${LIBSNARK_DIR}")
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # Common compilation flags and warning configuration
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Wextra -Wfatal-errors -pthread -maes -msse2 -msse3 -mpopcnt")

  if("${MULTICORE}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  endif()

  # Default optimizations flags (to override, use -DOPT_FLAGS=...)
  if("${OPT_FLAGS}" STREQUAL "")
    set(OPT_FLAGS "-ggdb3 -O2 -march=native -mtune=native")
  endif()
endif()

set(
  CMAKE_CXX_FLAGS
  "${CMAKE_CXX_FLAGS} ${OPT_FLAGS}"
)

include_directories(.)

add_subdirectory(depends)
add_subdirectory(lwe)
