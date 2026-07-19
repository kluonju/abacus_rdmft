# Unit tests for the basis-independent RDMFT optimisation core.  The core has
# no ABACUS dependencies, so the tested translation units are compiled directly
# into each test target (self-contained, no extra link libraries required).

AddTest(
    TARGET rdmft_core_xc_test
    SOURCES
        test_rdmft_xc.cpp
        ../rdmft_xc.cpp
)

AddTest(
    TARGET rdmft_core_occ_test
    SOURCES
        test_rdmft_occ.cpp
        ../rdmft_xc.cpp
        ../rdmft_occ_constraints.cpp
        ../rdmft_line_search.cpp
        ../rdmft_occ_optimizer.cpp
)

AddTest(
    TARGET rdmft_core_orbital_test
    SOURCES
        test_rdmft_orbital.cpp
        ../rdmft_occ_constraints.cpp
        ../rdmft_line_search.cpp
        ../rdmft_orbital_optimizer.cpp
)
