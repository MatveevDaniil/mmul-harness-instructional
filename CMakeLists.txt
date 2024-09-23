cmake_minimum_required(VERSION 3.14)
project(hw2 LANGUAGES C CXX)

set(CMAKE_BUILD_TYPE "Release")

find_package(BLAS REQUIRED)

# We require the GNU compiler for this assignment.
if (NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(Prg-Intel "PrgEnv-intel")
    set(Prg-Clang "PrgEnv-cray")
    message(WARNING
            "Must use GNU Compiler for submission. Make sure you ran:\n"
            "module swap ${Prg-${CMAKE_C_COMPILER_ID}} PrgEnv-gnu")
endif ()

# Common library target for benchmarking.
add_library(benchmark OBJECT benchmark.cpp)
target_compile_features(benchmark PRIVATE cxx_std_11)


set(BENCHMARKS_permutations
    "basic-ijk"
    "basic-ikj"
    "basic-jik"
    "basic-jki"
    "basic-kji"
    "basic-kij"
)
foreach(BENCHMARK ${BENCHMARKS_permutations})
    add_executable(benchmark-${BENCHMARK} dgemm-basic/dgemm-${BENCHMARK}.cpp)
    target_link_libraries(benchmark-${BENCHMARK} PRIVATE benchmark ${BLAS_LIBRARIES})
    target_compile_features(benchmark-${BENCHMARK} PRIVATE c_std_11 c_restrict)
    target_compile_options(benchmark-${BENCHMARK} PRIVATE -Wall -pedantic -march=native)
endforeach()


add_executable(benchmark-blas dgemm-blas.cpp)
target_link_libraries(benchmark-blas PRIVATE benchmark ${BLAS_LIBRARIES})
target_compile_features(benchmark-blas PRIVATE c_std_11 c_restrict)
target_compile_options(benchmark-blas PRIVATE -Wall -pedantic -march=native)


set(BENCHMARKS_blocked
    "blocked"
    "blocked-templated"
)
foreach(BENCHMARK ${BENCHMARKS_blocked})
    add_executable(benchmark-${BENCHMARK} dgemm-blocked/dgemm-${SOURCE_FILE}.cpp)
    add_library(bench-${BENCHMARK} OBJECT benchmark.cpp)
    target_compile_definitions(bench-${BENCHMARK} PRIVATE -DBLOCKED)
    target_compile_features(bench-${BENCHMARK} PRIVATE cxx_std_11)
    target_link_libraries(benchmark-${BENCHMARK} PRIVATE bench-${BENCHMARK} ${BLAS_LIBRARIES})
    target_compile_features(benchmark-${BENCHMARK} PRIVATE c_std_11 c_restrict)
    target_compile_options(benchmark-${BENCHMARK} PRIVATE -Wall -pedantic -march=native)
endforeach()


set(BENCHMARKS_strassen
    "strassen-1"
    "strassen-16"
    "strassen-64"
    "strassen-256"
)
foreach(BENCHMARK ${BENCHMARKS_strassen})
    add_executable(benchmark-${BENCHMARK} dgemm-strassen/dgemm-${BENCHMARK}.cpp)
    target_link_libraries(benchmark-${BENCHMARK} PRIVATE benchmark ${BLAS_LIBRARIES})
    target_compile_features(benchmark-${BENCHMARK} PRIVATE c_std_11 c_restrict)
    target_compile_options(benchmark-${BENCHMARK} PRIVATE -Wall -pedantic -march=native)
endforeach()
