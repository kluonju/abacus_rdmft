if(ENABLE_LCAO AND BUILD_TESTING)
    AddTest(
        TARGET rdmft_xc_test
        SOURCES test_rdmft_xc.cpp
    )
    AddTest(
        TARGET rdmft_occupation_test
        SOURCES test_rdmft_occupation.cpp
    )
    AddTest(
        TARGET rdmft_stiefel_test
        SOURCES test_rdmft_stiefel.cpp
    )
    AddTest(
        TARGET rdmft_optimizer_test
        SOURCES test_rdmft_optimizer.cpp
    )
endif()
