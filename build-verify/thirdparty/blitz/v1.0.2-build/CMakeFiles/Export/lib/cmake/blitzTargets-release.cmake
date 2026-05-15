#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "blitz" for configuration "Release"
set_property(TARGET blitz APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(blitz PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libblitz.so.1.0.2"
  IMPORTED_SONAME_RELEASE "libblitz.so.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS blitz )
list(APPEND _IMPORT_CHECK_FILES_FOR_blitz "${_IMPORT_PREFIX}/lib/libblitz.so.1.0.2" )

# Import target "blitz-static" for configuration "Release"
set_property(TARGET blitz-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(blitz-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libblitz.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS blitz-static )
list(APPEND _IMPORT_CHECK_FILES_FOR_blitz-static "${_IMPORT_PREFIX}/lib/libblitz.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
