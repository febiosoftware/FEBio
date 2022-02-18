cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
mark_as_advanced(FORCE CMAKE_INSTALL_PREFIX)

project(FEBio)
set(FEBIO_LIBS FEBioFluid FEBioLib FEBioMech FEBioMix FEBioOpt FEBioPlot FEBioTest FEBioXML FECore NumCore FEAMR FEBioRVE)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_FILES_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY AUTOGEN_TARGETS_FOLDER Autogen)

# Set a default build type if none was specified
set(default_build_type "Release")
 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

##### Set appropriate flag for operating system ####
if(WIN32)
	add_definitions(-DWIN32)
elseif(APPLE)
    add_definitions(-D__APPLE__)
else()
	add_definitions(-DLINUX)
endif()

##### Search for library and include directories #####

include(FindDependencies.cmake)

#### Check if paths are valid and find libs ####

function(findLib libName libDir libOut)
    # Optional arguments ARGV3-ARVG6 can be used as alernative names for the library
    if(WIN32)
        find_library(TEMP NAMES ${libName}.lib ${ARGV3}.lib ${ARGV4}.lib ${ARGV5}.lib ${ARGV6}.lib
            PATHS ${${libDir}} NO_DEFAULT_PATH)
    else()
        find_library(TEMP NAMES lib${libName}.a lib${ARGV3}.a lib${ARGV4}.a lib${ARGV5}.a lib${ARGV6}.a
            PATHS ${${libDir}} NO_DEFAULT_PATH)
    endif()
    
    if(TEMP)
        set(${libOut} ${TEMP} PARENT_SCOPE)
        unset(TEMP CACHE)
    else()
        if(WIN32)
            message(SEND_ERROR "Could not find ${libName}.lib. Check ${libName}.")
        elseif(APPLE)
            message(SEND_ERROR "Could not find lib${libName}.so, lib${libName}.a, or lib${libName}.dylib Check ${libDir}.")
        else()
            message(SEND_ERROR "Could not find lib${libName}.so, or lib${libName}.a. Check ${libDir}")
        endif()
        unset(TEMP CACHE)
    endif()
endfunction()

if(USE_MKL)
    set(MKL_LIB_NAMES mkl_intel_lp64 mkl_core mkl_intel_thread)
    foreach(name IN LISTS MKL_LIB_NAMES)
        findLib(${name} MKL_LIB_DIR MKL_LIB)
        list(APPEND MKL_LIBS ${MKL_LIB})
        unset(MMG_LIB)
    endforeach(name)
    
    if(NOT EXISTS ${MKL_OMP_LIB})
        message(SEND_ERROR "Could not find MKL OMP library. Check MKL_OMP_LIB.")
    endif()
endif()

if(USE_HYPRE)
	if(NOT EXISTS ${HYPRE_LIB})
		message(SEND_ERROR "Could not find HYPRE library. Check HYPRE_LIB.")
	endif()
	if(NOT EXISTS ${HYPRE_INC}/HYPRE_IJ_mv.h)
		message(SEND_ERROR "Could not find HYPRE_IJ_mv.h. Check HYPRE_INC.")
	endif()
endif()

if(USE_MMG)
	if(NOT EXISTS ${MMG_LIB})
		message(SEND_ERROR "Could not find MMG library. Check MMG_LIB.")
	endif()
	if(NOT EXISTS ${MMG_INC}/mmg/libmmg.h)
		message(SEND_ERROR "Could not find libmmg.h. Check MMG_INC.")
	endif()
endif()

if(USE_LEVMAR)
	if(NOT EXISTS ${LEVMAR_INC}/levmar.h)
		message(SEND_ERROR "Could not find levmar.h. Check LEVMAR_INC.")
	endif()
	if(NOT EXISTS ${LEVMAR_LIB})
		message(SEND_ERROR "Could not find Levmar library. Check LEVMAR_LIB.")
	endif()
endif()

if(USE_ZLIB)
	if(NOT EXISTS ${ZLIB_LIBRARY_RELEASE})
		message(SEND_ERROR "Could not find ZLIB. Check ZLIB_LIBRARY_RELEASE.")
	endif()
	if(NOT EXISTS ${ZLIB_INCLUDE_DIR}/zlib.h)
		message(SEND_ERROR "Could not find zlib.h. Check ZLIB_INCLUDE_DIR.")
	endif()
endif()

##### Set appropriate defines and includes #####
include_directories(${PROJECT_SOURCE_DIR})

if(USE_MKL)
	add_definitions(-DMKL_ISS)
	add_definitions(-DPARDISO)
	include_directories(${MKL_INC})
endif(USE_MKL)

if(USE_HYPRE)
	add_definitions(-DHYPRE)
	include_directories(${HYPRE_INC})
endif(USE_HYPRE)

if(USE_MMG)
	add_definitions(-DHAS_MMG)
	include_directories(${MMG_INC})
endif(USE_MMG)

if(USE_LEVMAR)
	add_definitions(-DHAVE_LEVMAR)
	include_directories(${LEVMAR_INC})
endif(USE_LEVMAR)

if(USE_ZLIB)
	add_definitions(-DHAVE_ZLIB)
	include_directories(${ZLIB_INCLUDE_DIR})
endif(USE_ZLIB)

##### Set up compiler flags #####
if(WIN32)
    add_compile_options(/MP /openmp)
elseif(APPLE)
	add_compile_options(-fPIC -w)
else()
	add_compile_options(-fPIC -fopenmp -w)
    
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++")
    
    # Extra compiler flags for intel compiler
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-intel")
    endif()
    
    set(CMAKE_BUILD_RPATH_USE_LINK_PATH FALSE)
    set(CMAKE_BUILD_RPATH $ORIGIN/../lib/)
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        add_compile_options(-no-intel-extensions)
endif()

# Add openMP flags for macOS if found.
if(APPLE)
    if(${OpenMP_C_FOUND})
        add_compile_options(-Xpreprocessor -fopenmp)
    endif()
endif()

##### Find Source Files #####

macro(findHdrSrc name)
	file(GLOB HDR_${name} "${name}/*.h")
	file(GLOB SRC_${name} "${name}/*.cpp")
endmacro()

findHdrSrc(FEBioStudio)

foreach(name IN LISTS FEBIO_LIBS)
	findHdrSrc(${name})
endforeach(name)

macro(addLib name)
    string(TOLOWER ${name} lname)
	add_library(${lname} ${HDR_${name}} ${SRC_${name}})
	set_property(TARGET ${lname} PROPERTY AUTOGEN_BUILD_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/AutoGen/${name}_autogen)
endmacro()

foreach(name IN LISTS FEBIO_LIBS)
	addLib(${name})
endforeach(name)

##### Set up executable compilation #####
file(GLOB SOURCES "FEBio3/*.cpp")
add_executable (febio3 ${SOURCES})

##### Set dev commit information #####

# Cross platform execute_process
macro(crossExecProcess OUT)
	if(WIN32)
		execute_process(COMMAND cmd.exe /c ${CMD}
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			RESULT_VARIABLE SUCCESS
			OUTPUT_VARIABLE ${OUT}
			OUTPUT_STRIP_TRAILING_WHITESPACE)
	else()
		execute_process(COMMAND ${CMD}
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			RESULT_VARIABLE SUCCESS
			OUTPUT_VARIABLE ${OUT}
			OUTPUT_STRIP_TRAILING_WHITESPACE)
	endif()
endmacro()

# Check for git and give warning if not available
if(EXISTS ${CMAKE_SOURCE_DIR}/.git)
    if(WIN32)
        set(CMD where git)
    else()
        set(CMD which git)
    endif()

        crossExecProcess(NULL)

    if(SUCCESS EQUAL 0)
        # Get branch name and check if we're on the develop branch
        set(CMD git rev-parse --abbrev-ref HEAD)
        crossExecProcess(BRANCH_NAME)

        string(COMPARE EQUAL ${BRANCH_NAME} "develop" ISDEVELOP)
        
        # Get the commit info and set the DEVCOMMIT macro
        if(SUCCESS EQUAL 0 AND ISDEVELOP)
            
            set(CMD git rev-parse --short=9 HEAD)
            crossExecProcess(DEVCOMMIT)
            
            if(${SUCCESS} EQUAL 0)
                target_compile_definitions(febiolib PUBLIC "DEVCOMMIT=\"${DEVCOMMIT}\"")
            endif()
        endif()
    else()
        message(WARNING "Could not find git in system PATH. Development version info will not be addded.")
    endif()
endif()

##### Linking options #####

# Link FEBio libraries
if(WIN32 OR APPLE)
	target_link_libraries(febio3 fecore febiolib febioplot febiomech 
		febiomix febioxml numcore febioopt febiotest febiofluid feamr febiorve)
else()
	target_link_libraries(febio3 -Wl,--start-group fecore febiolib febioplot febiomech 
		febiomix febioxml numcore febioopt febiotest febiofluid feamr febiorve -Wl,--end-group)
endif()

# Link LEVMAR
if(USE_LEVMAR)
	target_link_libraries(febio3 ${LEVMAR_LIB})
endif()

# Link HYPRE
if(USE_HYPRE)
	target_link_libraries(febio3 ${HYPRE_LIB})
endif()

# Link MKL
if(USE_MKL)
    if(WIN32 OR APPLE)
        target_link_libraries(febio3 ${MKL_LIBS} ${MKL_OMP_LIB})
    else()
        target_link_libraries(febio3 -Wl,--start-group ${MKL_LIBS} ${MKL_OMP_LIB} -Wl,--end-group)
    endif()
else()
    # If not using MKL, we still need OpenMP from the system.
    if(${OpenMP_C_FOUND})
        target_link_libraries(febio3 ${OpenMP_C_LIBRARIES})
    endif()
endif()

if(WIN32)
	target_link_libraries(febio3 psapi.lib ws2_32.lib)
else()
    target_link_libraries(febio3 -pthread -ldl)
endif()

# Link MMG
if(USE_MMG)
	target_link_libraries(febio3 ${MMG_LIB})
endif()

# Link ZLIB
if(USE_ZLIB)
	target_link_libraries(febio3 ${ZLIB_LIBRARY_RELEASE})
endif()

##### Create febio.xml #####
if(NOT EXISTS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/febio.xml)
    file(WRITE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/febio.xml "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>
<febio_config version=\"3.0\">
    <default_linear_solver type=\"pardiso\"></default_linear_solver>
</febio_config>")
endif()

if(USE_MKL)
    file(READ ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/febio.xml filedata)
    string(REGEX REPLACE "type=\"[a-z]*\"" "type=\"pardiso\"" filedata ${filedata})
    file(WRITE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/febio.xml ${filedata})
else()
    file(READ ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/febio.xml filedata)
    string(REGEX REPLACE "type=\"[a-z]*\"" "type=\"skyline\"" filedata "${filedata}")
    file(WRITE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/febio.xml "${filedata}")
endif()


get_target_property(LINKEDLIBS febio3 LINK_LIBRARIES)

