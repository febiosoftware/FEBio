# MKL - On Unix the compilervars.sh should be run to find the MKL libraries.
if(DEFINED ENV{MKLROOT})
    set(MKLROOT $ENV{MKLROOT} CACHE PATH "MKL root directory")
else()
    if(WIN32)
        set(MKLPATHS $ENV{ProgramFiles\(x86\)}/IntelSWTools $ENV{PROGRAMFILES}/Intel* $ENV{SystemDrive} $ENV{SystemDrive}/Intel*)
        set(MKLSUFFIXES "compilers_and_libraries/windows" "oneapi")
    elseif(APPLE)
        set(MKLPATHS /opt/intel /intel /usr/local/intel /usr/local/opt/intel)
        set(MKLSUFFIXES "compilers_and_libraries/mac" "oneapi")
    else()
        set(MKLPATHS /opt/intel /intel /usr/local/intel /usr/local/opt/intel $ENV{HOME}/intel $ENV{HOME}/*/intel)
        set(MKLSUFFIXES "compilers_and_libraries/linux" "oneapi")
    endif()
    
    find_file(MKLROOT mkl
		PATHS ${MKLPATHS}
		PATH_SUFFIXES ${MKLSUFFIXES}
		DOC "MKL root directory")
endif()

if(MKLROOT)
    if(${MKLROOT} MATCHES "oneapi")
        find_path(MKL_INC mkl.h 
            PATHS ${MKLROOT}/latest/include
            DOC "MKL include directory")
            
        find_library(MKL_CORE mkl_core
            PATHS ${MKLROOT}/latest/lib
            NO_DEFAULT_PATH)
            
        find_library(MKL_OMP_LIB 
            NAMES iomp5 iomp5md libiomp5md.lib
            PATHS ${MKLROOT}/../compiler/latest/*/compiler/lib
            NO_DEFAULT_PATH
            DOC "MKL OMP Library")
            
    else()
        find_path(MKL_INC mkl.h 
            PATHS ${MKLROOT}/include
            DOC "MKL include directory")
            
        find_library(MKL_CORE mkl_core
            PATHS ${MKLROOT}/lib
            PATH_SUFFIXES "intel64" "intel32"
            NO_DEFAULT_PATH)
            
        find_library(MKL_OMP_LIB 
            NAMES iomp5 iomp5md libiomp5md.lib
            PATHS ${MKLROOT}/lib ${MKLROOT}/../lib ${MKLROOT}/../compiler/lib
            PATH_SUFFIXES "intel64" "intel32"
            NO_DEFAULT_PATH
            DOC "MKL OMP Library")
    
    endif()
        
    if(MKL_CORE)
        get_filename_component(MKL_TEMP ${MKL_CORE} DIRECTORY)
        set(MKL_LIB_DIR ${MKL_TEMP} CACHE PATH "MKL Library directory")
        unset(MKL_TEMP)
        unset(MKL_CORE CACHE)
    else()
        set(MKL_LIB_DIR "MKL_LIB_DIR-NOTFOUND" CACHE PATH "MKL Library directory")
        unset(MKL_CORE CACHE)
    endif()
	
endif()

if(NOT WIN32)
    # OpenMP
    find_package(OpenMP QUIET)
endif()

if(MKL_INC AND MKL_LIB_DIR AND MKL_OMP_LIB)
	option(USE_MKL "Required for pardiso and iterative solvers" ON)
    mark_as_advanced(MKL_INC MKL_LIB_DIR MKL_OMP_LIB)
    set(OpenMP_C_FOUND true)
else()
	option(USE_MKL "Required for pardiso and iterative solvers" OFF)
endif()

# HYPRE
if(WIN32)
	find_path(HYPRE_INC HYPRE_IJ_mv.h
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
		PATH_SUFFIXES "include" "include/hypre" "src" "src/include" "src/hypre/include"
        DOC "HYPRE include directory")
	find_library(HYPRE_LIB HYPRE 
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
        PATH_SUFFIXES "src" "src/build" "src/mbuild" "/src/vs2017/Release"
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
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/* $ENV{HOMEPATH}/source/repos/*
		PATH_SUFFIXES "include" "include/mmg*" "src" "build" "build/include" "cmbuild/include"
        DOC "MMG include directory")
	find_library(MMG_LIB mmg3d 
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/* $ENV{HOMEPATH}/source/repos/*
        PATH_SUFFIXES "build/lib" "cmbuild/lib" "src/build/lib" "src/cmbuild/lib" "cmbuild/lib/Release"
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

# LEVMAR
if(WIN32)
	find_path(LEVMAR_INC levmar.h PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
		DOC "Levmar include directory")
	find_library(LEVMAR_LIB levmar PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
        PATH_SUFFIXES "vs2017/Release"
		DOC "Levmar library path")
else()
	find_path(LEVMAR_INC levmar.h PATHS /usr/local/ /opt/levmar* $ENV{HOME}/* $ENV{HOME}/*/*
		DOC "Levmar include directory")
	find_library(LEVMAR_LIB levmar PATHS /usr/local/ /opt/levmar* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "build" "cbuild" "cmbuild"
		DOC "Levmar library path")
endif()	

# SuperLU_MT
if (WIN32)
    find_path(SUPERLU_MT_INC slu_mt_ddefs.h PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
    DOC "SuperLU_MT include directory")
    find_library(SUPERLU_MT_LIB superlu PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
    DOC "SuperLU_MT library path")
else()
    find_path(SUPERLU_MT_INC slu_mt_ddefs.h PATHS /usr/local/include/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
    DOC "SuperLU_MT include directory")
    find_library(SUPERLU_MT_LIB superlu PATHS /usr/local/lib/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
    DOC "SuperLU_MT library path")
endif()

if(SUPERLU_MT_INC AND SUPERLU_MT_LIB)		
	option(USE_SUPERLU_MT "Option for using SuperLU_MT" ON)
    mark_as_advanced(SUPERLU_MT_INC SUPERLU_MT_LIB)
else()
	option(USE_SUPERLU_MT "Option for using SuperLU_MT" OFF)
    mark_as_advanced(CLEAR SUPERLU_MT_INC SUPERLU_MT_LIB)
endif()


if(LEVMAR_INC AND LEVMAR_LIB)		
	option(USE_LEVMAR "Required for optimization in FEBio" ON)
    mark_as_advanced(LEVMAR_INC LEVMAR_LIB)
else()
	option(USE_LEVMAR "Required for optimization in FEBio" OFF)
    mark_as_advanced(CLEAR LEVMAR_INC LEVMAR_LIB)
endif()

# PDL
if(WIN32)
	find_library(PDL_LIB libpardiso600-WIN-X86-64* 
        PATHS C::/Program\ Files/* $ENV{HOME}/* $ENV{HOME}/*/*
		DOC "PDL library path")
elseif(APPLE)
	find_library(PDL_LIB pardiso600-MACOS-X86-64 
        PATHS /usr/local* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "lib" "pardiso/lib" "pardiso-project/lib"
		DOC "PDL library path")
else()
	find_library(PDL_LIB pardiso600-GNU*-X86-64 
        PATHS /usr/local* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "lib" "pardiso/lib" "pardiso-project/lib"
		DOC "PDL library path")
endif()	

if(PDL_LIB)		
	option(USE_PDL "Required for pardiso-project use" ON)
    mark_as_advanced(PDL_LIB)
else()
	option(USE_PDL "Required for pardiso-project use" OFF)
    mark_as_advanced(CLEAR PDL_LIB)
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
