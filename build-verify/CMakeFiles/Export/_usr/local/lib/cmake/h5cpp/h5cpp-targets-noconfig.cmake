#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "libh5cpp-libdeflate-static" for configuration ""
set_property(TARGET libh5cpp-libdeflate-static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(libh5cpp-libdeflate-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/liblibh5cpp-libdeflate-static.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS libh5cpp-libdeflate-static )
list(APPEND _IMPORT_CHECK_FILES_FOR_libh5cpp-libdeflate-static "${_IMPORT_PREFIX}/lib/liblibh5cpp-libdeflate-static.a" )

# Import target "libh5cpp-zstd-static" for configuration ""
set_property(TARGET libh5cpp-zstd-static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(libh5cpp-zstd-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "ASM;C"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/liblibh5cpp-zstd-static.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS libh5cpp-zstd-static )
list(APPEND _IMPORT_CHECK_FILES_FOR_libh5cpp-zstd-static "${_IMPORT_PREFIX}/lib/liblibh5cpp-zstd-static.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
