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
#     FEBio::FEBioPlot
#     FEBio::FEBioTest
#     FEBio::FEBioXML
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------

set(_FEBIO_LIBS
    FECore
    FEBioMech
    FEBioMix
    FEBioFluid
    FEBioRVE
    FEBioPlot
    FEBioXML
    FEBioLib
    FEAMR
    FEBioOpt
    FEImgLib
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

    # Try to find the build directory by searching for fecore
    file(GLOB _MATCHED_DIRS LIST_DIRECTORIES true "${_FEBIO_PREFIX}/*build*")

    if(_MATCHED_DIRS)

        if(WIN32)
            set(_LIB_NAME "fecore.lib")
        elseif(APPLE)
            set(_LIB_NAME "libfecore.dylib")
        else()
            set(_LIB_NAME "libfecore.so")
        endif()


        foreach(dir IN LISTS _MATCHED_DIRS)
            find_file(_test
                NAMES ${_LIB_NAME}
                PATHS ${dir}
                PATH_SUFFIXES lib lib/Release Release/lib Debug/lib lib/Debug
                NO_DEFAULT_PATH            
            )

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

function(_febio_find_library out_release out_debug out_min_size out_rel_deb out_fallback libname)
    string(TOLOWER ${libname} _lib_lc)

    # Find release version
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

    # Find debug version
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

    # Find min size version
    find_library(_TEMP 
        NAMES ${_lib_lc}
        PATHS ${_FEBIO_LIB_DIR}
        PATH_SUFFIXES lib lib/MinSizeRel MinSizeRel/lib MinSizeRel
        NO_DEFAULT_PATH)

    string(REGEX MATCH "MinSizeRel" _is_min_size ${_TEMP})

    if(_is_min_size)
        set(${out_min_size} "${_TEMP}" PARENT_SCOPE)
    endif()
    unset(_TEMP CACHE)

    # Find rel with deb info version
    find_library(_TEMP 
        NAMES ${_lib_lc}
        PATHS ${_FEBIO_LIB_DIR}
        PATH_SUFFIXES lib lib/RelWithDebInfo RelWithDebInfo/lib RelWithDebInfo
        NO_DEFAULT_PATH)

    string(REGEX MATCH "RelWithDebInfo" _is_rel_deb ${_TEMP})

    if(_is_rel_deb)
        set(${out_rel_deb} "${_TEMP}" PARENT_SCOPE)
    endif()
    unset(_TEMP CACHE)

    # Fallback: find any version
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

foreach(lib IN LISTS _FEBIO_LIBS)
    if(TARGET "FEBio::${lib}")
        continue()
    endif()

    add_library("FEBio::${lib}" SHARED IMPORTED)

    set_target_properties("FEBio::${lib}" PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${_FEBIO_INCLUDE_DIR}"
    )

    _febio_find_library(_rel _dbg _min_size _rel_deb _fallback "${lib}")

    if(_rel)
        set_target_properties("FEBio::${lib}" PROPERTIES
            IMPORTED_LOCATION_RELEASE "${_rel}"
        )

        set_target_properties("FEBio::${lib}" PROPERTIES
            IMPORTED_IMPLIB_RELEASE "${_rel}"
        )

        # Use release configuration as default in case not all are found
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

        # Use debug configuration as default if release not found
        if(NOT _rel)
            set_target_properties("FEBio::${lib}" PROPERTIES
                IMPORTED_LOCATION "${_dbg}"
            )
            set_target_properties("FEBio::${lib}" PROPERTIES
                IMPORTED_IMPLIB "${_dbg}"
            )
        endif()
    endif()

    if(_min_size)
        set_target_properties("FEBio::${lib}" PROPERTIES
            IMPORTED_LOCATION_MINSIZEREL "${_min_size}"
        )
        set_target_properties("FEBio::${lib}" PROPERTIES
            IMPORTED_IMPLIB_MINSIZEREL "${_min_size}"
        )
    endif()

    if(_rel_deb)
        set_target_properties("FEBio::${lib}" PROPERTIES
            IMPORTED_LOCATION_RELWITHDEBINFO "${_rel_deb}"
        )
        set_target_properties("FEBio::${lib}" PROPERTIES
            IMPORTED_IMPLIB_RELWITHDEBINFO "${_rel_deb}"
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

unset(_FEBIO_LIBS)
unset(_FEBIO_PREFIX)
unset(_FEBIO_INCLUDE_DIR)
unset(_FEBIO_LIB_DIR)

set(FEBio_FOUND TRUE)