cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
mark_as_advanced(FORCE CMAKE_INSTALL_PREFIX)

project(FEBio)
set(FEBIO_LIBS FEBioFluid FEBioLib FEBioMech FEBioMix FEBioOpt FECore NumCore FEAMR FEBioRVE FEImgLib)
set(STATIC_FEBIO_LIBS FEBioPlot FEBioTest FEBioXML XML)

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
	add_definitions(-DWIN32 -DFECORE_DLL)
elseif(APPLE)
    add_definitions(-D__APPLE__)
    add_definitions(-DHAS_ACCEL)
else()
	add_definitions(-DLINUX)
    
    option(USE_STATIC_STDLIBS "Helps to create more portable binaries" ON)
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

if(USE_SUPERLU_MT)
	if(NOT EXISTS ${SUPERLU_MT_LIB})
		message(SEND_ERROR "Could not find SUPERLU_MT library. Check SUPERLU_MT_LIB.")
	endif()
	if(NOT EXISTS ${SUPERLU_MT_INC}/slu_mt_ddefs.h)
		message(SEND_ERROR "Could not find slu_mt_ddefs.h. Check SUPERLU_MT_INC.")
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

##### Set up compiler flags #####
if(WIN32)
    add_compile_options(/MP /wd4996 /wd4251 /wd4275)
elseif(APPLE)

else()
	add_compile_options(-fopenmp -w)
    
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

foreach(name IN LISTS STATIC_FEBIO_LIBS)
	findHdrSrc(${name})
endforeach(name)


macro(addLib name TYPE)
    string(TOLOWER ${name} lname)
	add_library(${lname} ${TYPE} ${HDR_${name}} ${SRC_${name}})
	set_property(TARGET ${lname} PROPERTY AUTOGEN_BUILD_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/AutoGen/${name}_autogen)
    
    if(NOT WIN32)
        set_property(TARGET ${lname} PROPERTY POSITION_INDEPENDENT_CODE ON)
    endif()
    
    if(APPLE)
        target_link_libraries(${lname} PRIVATE "-undefined dynamic_lookup")
    endif()
endmacro()

foreach(name IN LISTS FEBIO_LIBS)
	addLib(${name} SHARED)
endforeach(name)

foreach(name IN LISTS STATIC_FEBIO_LIBS)
	addLib(${name} STATIC)
endforeach(name)

##### Set up executable compilation #####
file(GLOB SOURCES "FEBio/*.cpp")
add_executable (febio4 ${SOURCES})

if(WIN32)
    target_compile_options(febio4 PRIVATE /openmp)
    target_compile_options(febiofluid PRIVATE /openmp)
    target_compile_options(febiolib PRIVATE /openmp)
    target_compile_options(febiomech PRIVATE /openmp)
    target_compile_options(febiomix PRIVATE /openmp)
    target_compile_options(febioopt PRIVATE /openmp)
    target_compile_options(fecore PRIVATE /openmp)
    target_compile_options(feamr PRIVATE /openmp)
    target_compile_options(febiorve PRIVATE /openmp)
    target_compile_options(feimglib PRIVATE /openmp)
    target_compile_options(febioplot PRIVATE /openmp)
    target_compile_options(febiotest PRIVATE /openmp)
    target_compile_options(febioxml PRIVATE /openmp)
    target_compile_options(xml PRIVATE /openmp)
endif()

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
        string(COMPARE EQUAL ${BRANCH_NAME} "febio4"  ISFEBIO4)
        
        # Get the commit info and set the DEVCOMMIT macro
        if(SUCCESS EQUAL 0 AND (ISDEVELOP OR ISFEBIO4))
            
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

##### Temp fix for Gerard to fix HYPRE's MPI include problem #####
set(EXTRA_INC "" CACHE STRING "Enter extra include paths separated by spaces")
set(EXTRA_LIB "" CACHE STRING "Enter extra library paths separated by spaces")
mark_as_advanced(EXTRA_INC)

##### Setup includes, defnitions, and linking options #####
include_directories(${PROJECT_SOURCE_DIR})

# Link LEVMAR
if(USE_LEVMAR)
	target_include_directories(febioopt PRIVATE ${LEVMAR_INC})
    target_compile_definitions(febioopt PRIVATE HAVE_LEVMAR)
    target_link_libraries(febioopt PRIVATE ${LEVMAR_LIB})
endif()

# Link HYPRE
if(USE_HYPRE)
    target_include_directories(numcore PRIVATE ${HYPRE_INC})
    target_compile_definitions(numcore PRIVATE HYPRE)
	target_link_libraries(numcore PRIVATE ${HYPRE_LIB})
endif()

# Link SuperLU_MT
if(USE_SUPERLU_MT)
	target_include_directories(numcore PRIVATE ${SUPERLU_MT_INC})
    target_compile_definitions(numcore PRIVATE HAVE_SUPERLU_MT)
    if(APPLE)
        target_compile_definitions(numcore PRIVATE _PTHREAD)
    endif()
    target_link_libraries(numcore PRIVATE ${SUPERLU_MT_LIB})
endif()

# Link MKL
if(USE_MKL)
    target_include_directories(numcore PRIVATE ${MKL_INC} ${OpenMP_C_LIBRARIES})
    target_compile_definitions(numcore PRIVATE MKL_ISS PARDISO)
    
    target_include_directories(febioopt PRIVATE ${MKL_INC})
	target_compile_definitions(febioopt PRIVATE MKL_ISS PARDISO)

    if(WIN32 OR APPLE)
        target_link_libraries(numcore PRIVATE ${MKL_OMP_LIB} ${MKL_LIBS})
        target_link_libraries(febioopt PRIVATE ${MKL_OMP_LIB} ${MKL_LIBS})
    else()
        target_link_libraries(numcore PRIVATE -Wl,--start-group ${MKL_OMP_LIB} ${MKL_LIBS} -Wl,--end-group)
        target_link_libraries(febioopt PRIVATE -Wl,--start-group ${MKL_OMP_LIB} ${MKL_LIBS} -Wl,--end-group)
    endif()
    
    if(FORCE_SYSTEM_OMP)
        target_link_libraries(numcore PRIVATE ${OpenMP_C_LIBRARIES})
        target_link_libraries(febioopt PRIVATE ${OpenMP_C_LIBRARIES})
    endif()
else()
    # If not using MKL, we still need OpenMP from the system.
    if(${OpenMP_C_FOUND})
        target_link_libraries(numcore PRIVATE ${OpenMP_C_LIBRARIES})
        target_link_libraries(febioopt PRIVATE ${OpenMP_C_LIBRARIES})
    endif()
endif()

# Link Pardiso-project
if(USE_PDL)
    target_compile_definitions(numcore PRIVATE PARDISODL)
    target_link_libraries(numcore PRIVATE ${PDL_LIB})

    # We still need OpenMP from the system.
    if(${OpenMP_C_FOUND})
        target_link_libraries(numcore PRIVATE ${OpenMP_C_LIBRARIES})
        target_link_libraries(febioopt PRIVATE ${OpenMP_C_LIBRARIES})
    endif()
endif()

# Link Accelerate framework
if (APPLE)
  target_link_libraries(numcore PRIVATE "-framework Accelerate")
endif()

# Link MMG
if(USE_MMG)
	target_include_directories(feamr PRIVATE ${MMG_INC})
    target_compile_definitions(feamr PRIVATE HAS_MMG)
    target_link_libraries(feamr PRIVATE ${MMG_LIB})
endif()

# Link ZLIB
if(USE_ZLIB)
    target_include_directories(febioplot PRIVATE ${ZLIB_INCLUDE_DIR})
    target_compile_definitions(febioplot PRIVATE HAVE_ZLIB)
	target_link_libraries(febioplot PRIVATE ${ZLIB_LIBRARY_RELEASE})
endif()

# Extra Includes
target_include_directories(febioopt PRIVATE ${EXTRA_INC})
target_include_directories(numcore PRIVATE ${EXTRA_INC})
target_include_directories(feamr PRIVATE ${EXTRA_INC})
target_include_directories(febiolib PRIVATE ${EXTRA_INC})

# Extra Libs
target_link_libraries(numcore PRIVATE ${EXTRA_LIB})


if(WIN32)
	target_link_libraries(febio4 psapi.lib ws2_32.lib)
else()
    target_link_libraries(febio4 -ldl)
endif()

# Link Libraries into FEBioLib
target_link_libraries(numcore PRIVATE fecore)
target_link_libraries(febioxml PRIVATE fecore)
target_link_libraries(febiotest PRIVATE fecore)
target_link_libraries(febiorve PRIVATE fecore febiomech febioxml febioplot xml)
target_link_libraries(febioplot PRIVATE fecore)
target_link_libraries(febioopt PRIVATE fecore febioxml xml)
target_link_libraries(febiomix PRIVATE febiomech fecore)
target_link_libraries(febiomech PRIVATE fecore)
target_link_libraries(febiofluid PRIVATE febiomix febiomech fecore)
target_link_libraries(feimglib PRIVATE fecore)
target_link_libraries(feamr PRIVATE fecore)

target_link_libraries(febiolib PRIVATE febioplot febiomech 
		febiomix febioxml numcore febioopt febiotest febiofluid feamr xml febiorve fecore feimglib)

if (WIN32)
    target_link_libraries(febiolib PRIVATE psapi.lib)
endif()

# Link FEBio libraries
if(WIN32 OR APPLE)
	target_link_libraries(febio4 fecore febiolib)
else()
    if(USE_STATIC_STDLIBS)
        target_link_libraries(numcore PRIVATE -static-libstdc++ -static-libgcc)
        target_link_libraries(febioxml PRIVATE -static-libstdc++ -static-libgcc)
        target_link_libraries(febiotest PRIVATE -static-libstdc++ -static-libgcc)
        target_link_libraries(febiorve PRIVATE -static-libstdc++ -static-libgcc)
        target_link_libraries(febioplot PRIVATE -static-libstdc++ -static-libgcc)
        target_link_libraries(febioopt PRIVATE -static-libstdc++ -static-libgcc)
        target_link_libraries(febiomix PRIVATE -static-libstdc++ -static-libgcc)
        target_link_libraries(febiomech PRIVATE -static-libstdc++ -static-libgcc)
        target_link_libraries(febiofluid PRIVATE -static-libstdc++ -static-libgcc)
        target_link_libraries(feimglib PRIVATE -static-libstdc++ -static-libgcc)
        target_link_libraries(feamr PRIVATE -static-libstdc++ -static-libgcc)
        target_link_libraries(febiolib PRIVATE -static-libstdc++ -static-libgcc)
        target_link_libraries(fecore PRIVATE -static-libstdc++ -static-libgcc)
        
        target_link_libraries(febio4 -static-libstdc++ -static-libgcc)
    endif()

	target_link_libraries(febio4 -Wl,--start-group fecore febiolib -Wl,--end-group)
    
    # Extra compiler flags for intel compiler
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        target_link_libraries(febio4 -static-intel)
    endif()
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

