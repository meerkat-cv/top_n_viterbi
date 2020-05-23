cmake_minimum_required(VERSION 3.8)

project(viterbi VERSION 0.0.1)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wfatal-errors -O2")

include_directories(
    "src/"
)

file(GLOB all_viterbi_src
    "src/*.cpp"
)


add_executable(viterbi_test ${all_viterbi_src})

target_link_libraries(viterbi_test)

add_library(viterbi SHARED ${all_viterbi_src})
