add_executable(LBMTest
    test.cpp
)

target_link_libraries(LBMTest PRIVATE
    gtest_main
    lbmLib
)

include(GoogleTest)
gtest_discover_tests(LBMTest)
