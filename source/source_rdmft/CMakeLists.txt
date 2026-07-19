# Basis-independent RDMFT optimisation core.
#
# This library contains no plane-wave / LCAO Hamiltonian code: it only depends
# on the C++ standard library, so it can be linked by both the plane-wave and
# LCAO RDMFT backends (see source_esolver/esolver_rdmft*).  It is compiled
# whenever ENABLE_RDMFT is set.

add_library(
    rdmft_core
    OBJECT
    rdmft_xc.cpp
    rdmft_occ_constraints.cpp
    rdmft_line_search.cpp
    rdmft_occ_optimizer.cpp
    rdmft_orbital_optimizer.cpp
    rdmft_driver.cpp
)

if(ENABLE_COVERAGE)
    add_coverage(rdmft_core)
endif()

if(BUILD_TESTING)
    add_subdirectory(test)
endif()
