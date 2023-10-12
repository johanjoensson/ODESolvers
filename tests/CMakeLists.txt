find_package(GTest REQUIRED)
# include_directories(${GTEST_INCLUDE_DIRS})

set(TEST_FILES
    # 	main_test.cpp
    numerov.cpp
    adams-moulton.cpp
)

add_executable(tests ${TEST_FILES})
target_link_libraries(tests GTest::gtest GTest::gtest_main ode-lib)
gtest_discover_tests(tests)
