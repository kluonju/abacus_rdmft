if(ENABLE_LCAO AND BUILD_TESTING)
    AddTest(
        TARGET rdmft_xc_test
        SOURCES test_rdmft_xc.cpp
    )
    AddTest(
        TARGET rdmft_occupation_test
        LIBS base ${math_libs} device parameter
        SOURCES test_rdmft_occupation.cpp
    )
    # oneAPI 2025.3 (IntelLLVM): icpx can ICE in HIR codegen on this TU (see rdmft/CMakeLists.txt).
    if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
        # -O2 -fno-vectorize still ICEs on oneAPI 2025.3; -O0 avoids HIR codegen crash.
        set_source_files_properties(
            test_rdmft_occupation.cpp
            PROPERTIES
            COMPILE_OPTIONS "-O0"
        )
    endif()
    AddTest(
        TARGET rdmft_optimizer_test
        SOURCES test_rdmft_optimizer.cpp
    )
    AddTest(
        TARGET rdmft_orbital_step_test
        LIBS ${math_libs}
        SOURCES test_rdmft_orbital_step.cpp
    )
    AddTest(
        TARGET rdmft_joint_test
        SOURCES test_rdmft_joint.cpp
    )
    AddTest(
        TARGET rdmft_occ_optimizer_test
        LIBS base ${math_libs} device parameter
        SOURCES test_rdmft_occ_optimizer.cpp
    )
endif()
