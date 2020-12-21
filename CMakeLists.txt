cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
mark_as_advanced(FORCE CMAKE_BUILD_TYPE CMAKE_INSTALL_PREFIX)

project(FEBio)

##### Set appropriate flag for operating system ####
if(WIN32)
	set(DEF -DWIN32)
elseif(APPLE)
    set(DEF -D__APPLE__)
else()
	set(DEF -DLINUX)
endif()

##### Check that paths are valid #####
if(USE_MKL)
	if(WIN32)
		if(NOT EXISTS ${MKL_LIB_DIR}/mkl_intel_lp64.lib)
			message(SEND_ERROR "Could not find mkl_intel_lp64.lib. Check MKL_LIB_DIR.")
		endif()
		if(NOT EXISTS ${MKL_LIB_DIR}/mkl_core.lib)
			message(SEND_ERROR "Could not find mkl_core.lib. Check MKL_LIB_DIR.")
		endif()
		if(NOT EXISTS ${MKL_LIB_DIR}/mkl_intel_thread.lib)
			message(SEND_ERROR "Could not find mkl_intel_thread.lib. Check MKL_LIB_DIR.")
		endif()
		if(NOT EXISTS ${OMP_LIB})
			message(SEND_ERROR "Could not find libiomp5md.lib. Check OMP_LIB.")
		endif()
		if(NOT EXISTS ${MKL_INC}/mkl.h)
			message(SEND_ERROR "Could not find mkl.h. Check MKL_INC.")
		endif()
	else()
		if(NOT EXISTS ${MKL_LIB_DIR}/libmkl_intel_lp64.a)
			message(SEND_ERROR "Could not find libmkl_intel_lp64.a. Check MKL_LIB_DIR.")
		endif()
		if(NOT EXISTS ${MKL_LIB_DIR}/libmkl_core.a)
			message(SEND_ERROR "Could not find libmkl_core.a. Check MKL_LIB_DIR.")
		endif()
		if(NOT EXISTS ${MKL_LIB_DIR}/libmkl_intel_thread.a)
			message(SEND_ERROR "Could not find libmkl_intel_thread.a. Check MKL_LIB_DIR.")
		endif()
		if(NOT EXISTS ${OMP_LIB})
			message(SEND_ERROR "Could not find libiomp5.so or libiomp5.dylib. Check OMP_LIB.")
		endif()
		if(NOT EXISTS ${MKL_INC}/mkl.h)
			message(SEND_ERROR "Could not find mkl.h. Check MKL_INC.")
		endif()
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

if(USE_TETGEN)
	if(NOT EXISTS ${TET_LIB})
		message(SEND_ERROR "Could not find TetGen library. Check TET_LIB.")
	endif()
	if(NOT EXISTS ${TET_INC}/tetgen.h)
		message(SEND_ERROR "Could not find tetgen.h. Check TET_INC.")
	endif()
endif()

if(USE_LEVMAR)
	if(NOT EXISTS ${LEVMAR_INC}/levmar.h)
		message(SEND_ERROR "Could not find levmar.h. Check LEVMAR_INC.")
	endif()
	if(NOT EXISTS ${LEVMAR_LIB})
		if(WIN32)
			message(SEND_ERROR "Could not find levmar.lib. Check LEVMAR_LIB.")
		else()
			message(SEND_ERROR "Could not find liblevmar.a. Check LEVMAR_LIB.")
		endif()
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


##### Search for libaries #####

# MKL - On Unix the compilervars.sh should be run to find the MKL libraries. 
if(WIN32)
	find_file(MKLROOT mkl
		PATHS $ENV{ProgramFiles\(x86\)}/Intel* $ENV{PROGRAMFILES}/Intel* $ENV{SystemDrive} $ENV{SystemDrive}/Intel*
		PATH_SUFFIXES "compilers_and_libraries/windows"
		DOC "MKL root directory")
		
	if(MKLROOT)
		find_path(MKL_INC mkl.h
			PATHS ${MKLROOT}/include
			DOC "MKL include directory")
		find_path(MKL_LIB_DIR mkl_core.lib 
			PATHS ${MKLROOT}/lib
			PATH_SUFFIXES "intel64"
			DOC "MKL Library directory")
		find_library(OMP_LIB iomp5md 
			#PATHS $ENV{ProgramFiles\(x86\)}/IntelSWTools/compilers_and_libraries/windows/compiler/lib/intel64
			PATHS ${MKLROOT}/../compiler/lib/ ${MKLROOT}/../lib/
			PATH_SUFFIXES "intel64" "intel32"
			DOC "OMP Library")		
	endif()
else()
	if(APPLE)
		find_file(MKLROOT mkl
			PATHS /opt/intel /intel /usr/local/intel /usr/local/opt/intel
			PATH_SUFFIXES "compilers_and_libraries/mac"
			DOC "MKL root directory")
	else()
		find_file(MKLROOT mkl
			PATHS /opt/intel /intel /usr/local/intel /usr/local/opt/intel
			PATH_SUFFIXES "compilers_and_libraries/linux"
			DOC "MKL root directory")
	endif()

	if(MKLROOT)
		find_path(MKL_INC mkl.h
			PATHS ${MKLROOT}/include
			DOC "MKL include directory")
		find_path(MKL_LIB_DIR libmkl_core.a 
			PATHS ${MKLROOT}/lib
			PATH_SUFFIXES "intel64"
			DOC "MKL Library directory")
		find_library(OMP_LIB iomp5
			PATHS ${MKLROOT}/../lib
			PATH_SUFFIXES "intel64" "intel32"
			NO_DEFAULT_PATH
			DOC "OMP Library")		
	endif()
	
endif()

if(MKL_INC AND MKL_LIB_DIR AND OMP_LIB)
	option(USE_MKL "Required for pardiso and iterative solvers" ON)
else()
	option(USE_MKL "Required for pardiso and iterative solvers" OFF)
endif()

# HYPRE
if(WIN32)
	find_path(HYPRE_INC HYPRE.h 
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
		PATH_SUFFIXES "include" "include/hypre" "src" "src/include" "src/hypre/include"
        DOC "HYPRE include directory")
	find_library(HYPRE_LIB HYPRE 
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
        PATH_SUFFIXES "src" "src/build" "src/mbuild"
		DOC "HYPRE library path")
else()
	find_path(HYPRE_INC HYPRE_IJ_mv.h
        PATHS /opt/hypre* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "include" "include/hypre" "src" "src/include" "src/hypre/include"
		DOC "HYPRE include directory")
	find_library(HYPRE_LIB HYPRE 
        PATHS /opt/hypre* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "src" "src/build" "src/cbuild"
		DOC "HYPRE library path")
endif()	

if(HYPRE_INC AND HYPRE_LIB)		
	option(USE_HYPRE "Required for HYPRE solver" ON)
    mark_as_advanced(HYPRE_INC HYPRE_LIB)
else()
	option(USE_HYPRE "Required for HYPRE solver" OFF)
    mark_as_advanced(CLEAR HYPRE_INC HYPRE_LIB)
endif()

# MMG
if(WIN32)
	find_path(MMG_INC mmg/mmg3d/libmmg3d.h
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
		PATH_SUFFIXES "include" "include/mmg*" "src" "build" "build/include"
        DOC "MMG include directory")
	find_library(MMG_LIB mmg3d 
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
        PATH_SUFFIXES "build/lib" "cmbuild/lib" "src/build/lib" "src/cmbuild/lib"
		DOC "MMG library path")
else()
	find_path(MMG_INC mmg/mmg3d/libmmg3d.h
        PATHS /opt/hypre* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "include" "include/mmg" "build" "build/include" "cbuild" "cbuild/include" "src" 
		DOC "MMG include directory")
	find_library(MMG_LIB mmg3d 
        PATHS /opt/mmg* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "build/lib" "cbuild/lib" "src/build/lib" "src/cbuild/lib"
		DOC "MMG library path")
endif()	

if(MMG_INC AND MMG_LIB)		
	option(USE_MMG "Required for MMG use" ON)
    mark_as_advanced(MMG_INC MMG_LIB)
else()
	option(USE_MMG "Required for MMG use" OFF)
    mark_as_advanced(CLEAR MMG_INC MMG_LIB)
endif()

# TETGEN
if(WIN32)
	find_path(TET_INC tetgen.h
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
		PATH_SUFFIXES "include" "include/tetgen*" "src" "build" "build/include"
        DOC "TetGen include directory")
	find_library(TET_LIB tet 
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
        PATH_SUFFIXES "build/lib" "cmbuild/lib" "src/build/lib" "src/cmbuild/lib"
		DOC "TetGen library path")
else()
	find_path(TET_INC tetgen.h
        PATHS /opt/tetgen* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "include" "include/tetgen*" "src" "build" "build/include"
		DOC "TetGen include directory")
	find_library(TET_LIB tet 
        PATHS /opt/tetgen* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "build/lib" "cmbuild/lib" "src/build/lib" "src/cmbuild/lib"
		DOC "TetGen library path")
endif()	

if(TET_INC AND TET_LIB)		
	option(USE_TETGEN "Required for adaptive remeshing" ON)
    mark_as_advanced(TET_INC TET_LIB)
else()
	option(USE_TETGEN "Required for adaptive remeshing" OFF)
    mark_as_advanced(CLEAR TET_INC TET_LIB)
endif()

# LEVMAR
if(WIN32)
	find_path(LEVMAR_INC levmar.h PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
		DOC "Levmar include directory")
	find_library(LEVMAR_LIB levmar PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
		DOC "Levmar library path")
else()
	find_path(LEVMAR_INC levmar.h PATHS /usr/local/ /opt/levmar* $ENV{HOME}/* $ENV{HOME}/*/*
		DOC "Levmar include directory")
	find_library(LEVMAR_LIB levmar PATHS /usr/local/ /opt/levmar* $ENV{HOME}/* $ENV{HOME}/*/*
		DOC "Levmar library path")
endif()	

if(LEVMAR_INC AND LEVMAR_LIB)		
	option(USE_LEVMAR "Required for optimization in FEBio" ON)
    mark_as_advanced(LEVMAR_INC LEVMAR_LIB)
else()
	option(USE_LEVMAR "Required for optimization in FEBio" OFF)
    mark_as_advanced(CLEAR LEVMAR_INC LEVMAR_LIB)
endif()

# ZLIB
include(FindZLIB)

if(ZLIB_INCLUDE_DIR AND ZLIB_LIBRARY_RELEASE)		
	option(USE_ZLIB "Required for compressing xplt files" ON)
    mark_as_advanced(ZLIB_INCLUDE_DIR ZLIB_LIBRARY_RELEASE)
else()
	option(USE_ZLIB "Required for compressing xplt files" OFF)
    mark_as_advanced(CLEAR ZLIB_INCLUDE_DIR ZLIB_LIBRARY_RELEASE)
endif()

##### Set appropriate defines and includes #####
set(INC ${PROJECT_SOURCE_DIR})

if(USE_MKL)
	set(DEF ${DEF} -DMKL_ISS)
	set(DEF ${DEF} -DPARDISO)
	set(INC ${INC} ${MKL_INC})
endif(USE_MKL)

if(USE_HYPRE)
	set(DEF ${DEF} -DHYPRE)
	set(INC ${INC} ${HYPRE_INC})
endif(USE_HYPRE)

if(USE_MMG)
	set(DEF ${DEF} -DHAS_MMG)
	set(INC ${INC} ${MMG_INC})
endif(USE_MMG)

if(USE_TETGEN)
	set(DEF ${DEF} -DTETLIBRARY)
	set(INC ${INC} ${TET_INC})
endif(USE_TETGEN)

if(USE_LEVMAR)
	set(DEF ${DEF} -DHAVE_LEVMAR)
	set(INC ${INC} ${LEVMAR_INC})
endif(USE_LEVMAR)

if(USE_ZLIB)
	set(DEF ${DEF} -DHAVE_ZLIB)
	set(INC ${INC} ${ZLIB_INCLUDE_DIR})
endif(USE_ZLIB)

set(DEF ${DEF} -DNDEBUG)

add_definitions(${DEF})
include_directories(${INC})

##### Set up compiler flags #####
if(WIN32)
elseif(APPLE)
	add_definitions(-O3 -fPIC -std=c++11 -static-libstdc++ -w)
else()
	add_definitions(-O3 -fPIC -fopenmp -std=c++11 -static-libstdc++ -w)
    
    set(CMAKE_BUILD_RPATH_USE_LINK_PATH FALSE)
    set(CMAKE_BUILD_RPATH $ORIGIN/../lib/)
endif()

# Extra compiler flags for intel compiler
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	add_definitions(-static-intel -no-intel-extensions)
endif()

# Properly use OMP if llvm has been installed from brew
if(APPLE AND "${CMAKE_CXX_COMPILER}" MATCHES "opt")
	add_definitions(-fopenmp)
endif()

##### Include FEBio libraries #####
add_subdirectory (FEBioFluid)
add_subdirectory (FEBioLib)
add_subdirectory (FEBioMech)
add_subdirectory (FEBioMix)
add_subdirectory (FEBioOpt)
add_subdirectory (FEBioPlot)
add_subdirectory (FEBioTest)
add_subdirectory (FEBioXML)
add_subdirectory (FECore)
add_subdirectory (NumCore)

##### Set up executable compilation #####
file(GLOB SOURCES "FEBio3/*.cpp")
add_executable (febio3 ${SOURCES})

##### Set dev commit information #####

execute_process(COMMAND git rev-parse --abbrev-ref HEAD
    RESULT_VARIABLE SUCCESS
    OUTPUT_VARIABLE BRANCH_NAME)

string(STRIP ${BRANCH_NAME} BRANCH_NAME)
string(COMPARE EQUAL ${BRANCH_NAME} "develop" ISDEVELOP)
    
if(SUCCESS EQUAL 0 AND ISDEVELOP)
    execute_process(COMMAND git rev-parse --short HEAD
        RESULT_VARIABLE SUCCESS
        OUTPUT_VARIABLE DEVCOMMIT)
    
    string(STRIP ${DEVCOMMIT} DEVCOMMIT)
    
    if(${SUCCESS} EQUAL 0)
        target_compile_definitions(febiolib PUBLIC "DEVCOMMIT=\"${DEVCOMMIT}\"")
    endif()
endif()


##### Linking options #####

# Link FEBio libraries
if(WIN32 OR APPLE)
	target_link_libraries(febio3 fecore febiolib febioplot febiomech 
		febiomix febioxml numcore febioopt febiotest febiofluid)
else()
	target_link_libraries(febio3 -Wl,--start-group fecore febiolib febioplot febiomech 
		febiomix febioxml numcore febioopt febiotest febiofluid -Wl,--end-group)
endif()

# Link LEVMAR
if(USE_LEVMAR)
	target_link_libraries(febio3 ${LEVMAR_LIB})
endif(USE_LEVMAR)

# Link HYPRE
if(USE_HYPRE)
	target_link_libraries(febio3 ${HYPRE_LIB})
endif(USE_HYPRE)

# Link MKL
if(USE_MKL)
	if(WIN32)
		target_link_libraries(febio3 ${MKL_LIB_DIR}/mkl_intel_lp64.lib ${MKL_LIB_DIR}/mkl_core.lib
			${MKL_LIB_DIR}/mkl_intel_thread.lib)
		target_link_libraries(febio3 ${OMP_LIB})
	elseif(APPLE)
		target_link_libraries(febio3 ${MKL_LIB_DIR}/libmkl_intel_lp64.a ${MKL_LIB_DIR}/libmkl_core.a
			${MKL_LIB_DIR}/libmkl_intel_thread.a)
		target_link_libraries(febio3 ${OMP_LIB} -pthread -ldl)
    else()
        target_link_libraries(febio3 -Wl,--start-group ${MKL_LIB_DIR}/libmkl_intel_lp64.a ${MKL_LIB_DIR}/libmkl_core.a
			${MKL_LIB_DIR}/libmkl_intel_thread.a -Wl,--end-group)
		target_link_libraries(febio3 ${OMP_LIB} -pthread -ldl)
	endif()
endif()

# Link MMG
if(USE_MMG)
	target_link_libraries(febio3 ${MMG_LIB})
endif(USE_MMG)

# Link TETGEN
if(USE_TETGEN)
	target_link_libraries(febio3 ${TET_LIB})
endif(USE_TETGEN)

# Link ZLIB
if(USE_ZLIB)
	target_link_libraries(febio3 ${ZLIB_LIBRARY_RELEASE})
endif(USE_ZLIB)

##### Create febio.xml #####
if(USE_MKL)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/febio.xml "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>
<febio_config version=\"3.0\">
	<default_linear_solver type=\"pardiso\"></default_linear_solver>
</febio_config>
")
else()
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/febio.xml "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>
<febio_config version=\"3.0\">
	<default_linear_solver type=\"skyline\"></default_linear_solver>
</febio_config>
")
endif()


