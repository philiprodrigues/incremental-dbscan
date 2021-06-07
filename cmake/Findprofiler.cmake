# Attempt to find perftools' libprofiler via pkgconfig

find_package(PkgConfig)
pkg_check_modules(PC_profiler QUIET REQUIRED libprofiler)

find_path(profiler_INCLUDE_DIR
  NAMES gperftools/profiler.h
  PATHS ${PC_profiler_INCLUDE_DIRS}
  PATH_SUFFIXES gperftools
)

find_library(profiler_LIBRARY
  NAMES profiler
  PATHS ${PC_profiler_LIBRARY_DIRS}
  )

set(profiler_VERSION ${PC_profiler_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(profiler
  FOUND_VAR profiler_FOUND
  REQUIRED_VARS
    profiler_LIBRARY
    profiler_INCLUDE_DIR
  VERSION_VAR profiler_VERSION
  )

if(profiler_FOUND AND NOT TARGET profiler::profiler)
  add_library(profiler::profiler UNKNOWN IMPORTED)
  set_target_properties(profiler::profiler PROPERTIES
    IMPORTED_LOCATION "${profiler_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_profiler_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${profiler_INCLUDE_DIR}"
  )
endif()
