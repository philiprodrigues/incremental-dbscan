cmake_minimum_required(VERSION 3.17)
project(incremental-dbscan VERSION 0.0.1)

set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_MODULE_PATH})

find_package(ROOT 6.22 CONFIG REQUIRED)
find_package(profiler MODULE)
if(profiler_FOUND)
  add_compile_definitions(HAVE_PROFILER)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native")

add_executable(run_dbscan run_dbscan.cxx Hit.cpp draw_clusters.cpp dbscan_orig.cpp dbscan.cpp)
target_link_libraries(run_dbscan PUBLIC ROOT::Core ROOT::Graf ROOT::Rint ROOT::Gpad)
if(profiler_FOUND)
  target_link_libraries(run_dbscan PUBLIC profiler::profiler)
endif()
