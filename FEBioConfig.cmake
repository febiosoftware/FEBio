cmake_minimum_required(VERSION 3.15)

# ------------------------------------------------------------------------------
# FEBioConfig.cmake
#
# Exposes one imported target per FEBio library.
#
# Targets:
#     FEBio::FEAMR
#     FEBio::FEBioLib
#     FEBio::FEBioMech
#     FEBio::FEBioOpt
#     FEBio::FEBioRVE
#     FEBio::FECore
#     FEBio::FEImgLib
#     FEBio::NumCore
#
#     FEBio::FEBioPlot     (static)
#     FEBio::FEBioTest     (static)
#     FEBio::FEBioXML        (static)
#     FEBio::XML                 (static)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------

set(_FEBIO_SHARED_LIBS
    FEAMR
    FEBioLib
    FEBioMech
    FEBioOpt
    FEBioRVE
    FECore
    FEImgLib
    NumCore
)

set(_FEBIO_STATIC_LIBS
    FEBioPlot
    FEBioTest
    FEBioXML
    XML
)

# ------------------------------------------------------------------------------
# Detect in-tree build
# ------------------------------------------------------------------------------

get_filename_component(_FEBIO_CONFIG_DIR
    "${CMAKE_CURRENT_LIST_FILE}" PATH
)

set(_FEBIO_IS_SOURCE_TREE TRUE)
if(EXISTS "${_FEBIO_CONFIG_DIR}/../../../include")
        set(_FEBIO_IS_SOURCE_TREE FALSE)
endif()

if(_FEBIO_IS_SOURCE_TREE)
    # Config file is at the root of the FEBio source tree
    set(_FEBIO_PREFIX "${_FEBIO_CONFIG_DIR}")
    set(_FEBIO_INCLUDE_DIR "${_FEBIO_PREFIX}")

    # This is the default build output directory
    set(_FEBIO_LIB_DIR "${_FEBIO_PREFIX}/build")

    # Try to find the build directory by searching for febioxml
    file(GLOB _MATCHED_DIRS LIST_DIRECTORIES true "${_FEBIO_PREFIX}/*build*")

    if(_MATCHED_DIRS)
        foreach(dir IN LISTS _MATCHED_DIRS)
            find_library(_test 
                NAMES febioxml
                PATHS ${dir}
                PATH_SUFFIXES lib lib/Release Release/lib Debug/lib lib/Debug
                NO_DEFAULT_PATH)

            if(_test)
                set(_FEBIO_LIB_DIR "${dir}")
                break()
            endif()
        endforeach()

        unset(_test CACHE)
    endif()
else()
    # Installed SDK: <prefix>/lib/cmake/FEBio/FEBioConfig.cmake
    get_filename_component(_FEBIO_PREFIX
        "${_FEBIO_CONFIG_DIR}/../../.." REALPATH
    )
    set(_FEBIO_INCLUDE_DIR "${_FEBIO_PREFIX}/include")
    set(_FEBIO_LIB_DIR "${_FEBIO_PREFIX}/lib")
endif()

# ------------------------------------------------------------------------------
# Helper to locate libraries
# ------------------------------------------------------------------------------

function(_febio_find_library out_release out_debug out_fallback libname is_static)
    string(TOLOWER ${libname} _lib_lc)

    find_library(_TEMP
        NAMES ${_lib_lc}
        PATHS ${_FEBIO_LIB_DIR}
        PATH_SUFFIXES lib lib/Release Release/lib Release
        NO_DEFAULT_PATH)

    string(REGEX MATCH "Release" _is_release ${_TEMP})

    if(_is_release)
        set(${out_release} ${_TEMP} PARENT_SCOPE)
    endif()
    unset(_TEMP CACHE)

    find_library(_TEMP 
        NAMES ${_lib_lc}
        PATHS ${_FEBIO_LIB_DIR}
        PATH_SUFFIXES lib lib/Debug Debug/lib Debug
        NO_DEFAULT_PATH)

    string(REGEX MATCH "Debug" _is_debug ${_TEMP})

    if(_is_debug)
        set(${out_debug} "${_TEMP}" PARENT_SCOPE)
    endif()
    unset(_TEMP CACHE)

    find_library(_TEMP 
        NAMES ${_lib_lc}
        PATHS ${_FEBIO_LIB_DIR}
        PATH_SUFFIXES lib
        NO_DEFAULT_PATH)

    set(${out_fallback} "${_TEMP}" PARENT_SCOPE)
    unset(_TEMP CACHE)
endfunction()

# ------------------------------------------------------------------------------
# Create imported targets
# ------------------------------------------------------------------------------

foreach(lib IN LISTS _FEBIO_SHARED_LIBS _FEBIO_STATIC_LIBS)
    if(TARGET "FEBio::${lib}")
        continue()
    endif()

    list(FIND _FEBIO_STATIC_LIBS "${lib}" _static_index)
    if(_static_index GREATER -1)
        set(_lib_type STATIC)
        set(_is_static TRUE)
    else()
        set(_lib_type SHARED)
        set(_is_static FALSE)
    endif()

    add_library("FEBio::${lib}" ${_lib_type} IMPORTED)

    set_target_properties("FEBio::${lib}" PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${_FEBIO_INCLUDE_DIR}"
    )

    _febio_find_library(_rel _dbg _fallback "${lib}" "${_is_static}")

    if(_rel)
        set_target_properties("FEBio::${lib}" PROPERTIES
            IMPORTED_LOCATION_RELEASE "${_rel}"
        )

        set_target_properties("FEBio::${lib}" PROPERTIES
            IMPORTED_IMPLIB_RELEASE "${_rel}"
        )

        set_target_properties("FEBio::${lib}" PROPERTIES
            IMPORTED_LOCATION "${_rel}"
        )

        set_target_properties("FEBio::${lib}" PROPERTIES
            IMPORTED_IMPLIB "${_rel}"
        )        
    endif()
   
    if(_dbg)
        set_target_properties("FEBio::${lib}" PROPERTIES
            IMPORTED_LOCATION_DEBUG "${_dbg}"
        )
        set_target_properties("FEBio::${lib}" PROPERTIES
            IMPORTED_IMPLIB_DEBUG "${_dbg}"
        )
    endif()

    if(_fallback)
        set_target_properties("FEBio::${lib}" PROPERTIES
            IMPORTED_LOCATION "${_fallback}"
        )
        set_target_properties("FEBio::${lib}" PROPERTIES
            IMPORTED_IMPLIB "${_fallback}"
        )
    endif()

    if(NOT _rel AND NOT _dbg AND NOT _fallback)
        message(
            "FEBio library '${lib}' could not be found.\n"
            "Expected under: ${_FEBIO_PREFIX}"
        )

        set(FEBio_FOUND FALSE)
        return()
    endif()

    unset(_rel)
    unset(_dbg)
    unset(_fallback)

endforeach()

unset(_FEBIO_SHARED_LIBS)
unset(_FEBIO_STATIC_LIBS)
unset(_FEBIO_PREFIX)
unset(_FEBIO_INCLUDE_DIR)

set(FEBio_FOUND TRUE)