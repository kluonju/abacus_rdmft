if(ENABLE_LCAO AND BUILD_TESTING)
    add_executable(test_rdmft_xc test_rdmft_xc.cpp)
    target_link_libraries(test_rdmft_xc PRIVATE GTest::gtest_main)
    add_test(NAME test_rdmft_xc COMMAND test_rdmft_xc)

    add_executable(test_rdmft_occupation test_rdmft_occupation.cpp)
    target_link_libraries(test_rdmft_occupation PRIVATE GTest::gtest_main)
    add_test(NAME test_rdmft_occupation COMMAND test_rdmft_occupation)

    add_executable(test_rdmft_stiefel test_rdmft_stiefel.cpp)
    target_link_libraries(test_rdmft_stiefel PRIVATE GTest::gtest_main)
    add_test(NAME test_rdmft_stiefel COMMAND test_rdmft_stiefel)

    add_executable(test_rdmft_optimizer test_rdmft_optimizer.cpp)
    target_link_libraries(test_rdmft_optimizer PRIVATE GTest::gtest_main)
    add_test(NAME test_rdmft_optimizer COMMAND test_rdmft_optimizer)
endif()
