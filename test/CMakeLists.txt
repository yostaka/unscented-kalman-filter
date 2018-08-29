project(unitTests)

set(sources src/ukf.cpp src/tools.cpp)

add_executable(unitTests
        Tests.cpp)

target_link_libraries(unitTests
        gtest
        gtest_main
        pthread)
