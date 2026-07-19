list(APPEND objects
    esolver.cpp
    esolver_ks.cpp
    esolver_fp.cpp
    esolver_ks_pw.cpp
    esolver_ks_lcaopw.cpp
    esolver_sdft_pw.cpp
    esolver_lj.cpp
    esolver_dp.cpp
    esolver_nep.cpp
    esolver_of.cpp
    esolver_of_tddft.cpp
    esolver_of_interface.cpp
    esolver_of_tool.cpp
    pw_others.cpp
    esolver_dfpt_pw.cpp
)
if(ENABLE_LCAO)
  list(APPEND objects
      esolver_ks_lcao.cpp
      esolver_ks_lcao_tddft.cpp
      esolver_gets.cpp
      lcao_others.cpp
      esolver_dm2rho.cpp
      esolver_double_xc.cpp
  )
  if(ENABLE_RDMFT)
    # New modular RDMFT esolver (LCAO/LibRI backend around rdmft::EnergyGradient).
    list(APPEND objects esolver_rdmft.cpp)
  endif()
endif()

# Plane-wave RDMFT esolver (ACE exchange backend) is built when RDMFT is on
# without LCAO, i.e. the ACE backend replaces the LibRI backend.
if(ENABLE_RDMFT AND NOT ENABLE_LCAO)
  list(APPEND objects esolver_rdmft_pw.cpp)
endif()

add_library(
    esolver
    OBJECT
    ${objects}
        ../source_pw/module_pwdft/exx_helper.cpp
        ../source_pw/module_pwdft/exx_helper.h
)

if(ENABLE_RDMFT)
  target_compile_definitions(esolver PRIVATE __RDMFT)
endif()

if(ENABLE_COVERAGE)
  add_coverage(esolver)
endif()

if(BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()

