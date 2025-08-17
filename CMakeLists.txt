cmake_minimum_required(VERSION 3.20)

project(MMFT-LbmCore
    DESCRIPTION "A lightweight LBM solver that can handle CFD simulators in MMFT-Simulator"
    HOMEPAGE_URL "https://github.com/cda-tum/mmft-LbmCore.git"
    LANGUAGES CXX C
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# download external libraries
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.17.0
)
FetchContent_MakeAvailable(googletest)

# add library
set(TARGET_NAME lbmLib)
add_library(${TARGET_NAME})

# add sources
add_subdirectory(src)

# main executable
add_executable(MMFTLBM)
target_link_libraries(MMFTLBM PUBLIC ${TARGET_NAME})
target_sources(MMFTLBM PUBLIC src/main.cpp)

# create tests
option(TEST "Enable unit tests" OFF)
if(TEST)
  include(FetchContent)
  FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG v1.17.0
  )
  FetchContent_MakeAvailable(googletest)

  enable_testing()
  add_subdirectory(tests)
endif()
