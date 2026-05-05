if(ENABLE_LCAO)
    add_library(
        rdmft
        OBJECT
        rdmft.cpp
        rdmft_tools.cpp
        rdmft_pot.cpp
        update_state_rdmft.cpp
        rdmft_energy_gradient.cpp
        rdmft_restart_io.cpp
        rdmft_solver.cpp
    )
endif()


    # oneAPI 2025.3 can ICE in HIR codegen on RDMFTSolver::solve at higher -O levels.
    # Keep optimization high globally, but compile this TU conservatively as a workaround.
    if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
        set_source_files_properties(
            rdmft_solver.cpp
            PROPERTIES
            COMPILE_OPTIONS "-O2;-fno-vectorize"
        )
    endif()

if (BUILD_TESTING)
    add_subdirectory(test)
endif()
