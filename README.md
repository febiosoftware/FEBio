<img src="Documentation/febio_logo_full.png" href="https://gibboncode.org" alt="FEBio" width="75%">
<!-- ![febio logo](Documentation/febio_logo_full.png) -->

FEBio is a nonlinear finite element (FE) solver that is specifically designed for biomechanical applications. It offers modeling scenarios, constitutive models and boundary conditions that are relevant to many research areas in biomechanics. All features can be used together seamlessly, giving the user a powerful tool for solving 3D problems in computational biomechanics. The software is open-source, and pre-compiled executables for Windows, OS-X and Linux platforms are available.

Executables for FEBio can be downloaded from https://febio.org/febio/febio-downloads/.  Please inform us of publications that use FEBio in research.  Information can be found on the Publications tab.  

Support forums can be found at https://forums.febio.org/.

[![DOI](https://img.shields.io/badge/Citation-DOI:10.1115/1.4005694-green.svg)](https://dx.doi.org/10.1115%2F1.4005694)

[![License](https://img.shields.io/badge/License-MIT-orange.svg)](LICENSE)

### Table of contents
- [FEBio BUILD GUIDE](#Build)  
- [Contributing](#Contributing)  

# FEBio BUILD GUIDE <a name="Build"></a>

## Third Party Libraries

FEBio relies on several third party libraries. It is possible to compile FEBio without these libraries but it will be lacking some features.

* FEBio requires Intel MKL in order to utilize the Pardiso linear solver. This library can be downloaded from Intel's website: https://software.intel.com/en-us/mkl/. In order to compile using this library the MKL_ISS compiler flag must be set. Alternatively, FEBio's source comes with the default linear solver Skyline, which will be used in the absence of Intel MKL. However, the Pardiso solver is significantly faster and more memory-efficient than the Skyline solver, and it is strongly recommended that the Pardiso solver be used.  

* The Lourakis levmar routine is required by FEBio in order for it to perform its parameter optimization functions. The source for this library can be downloaded from http://users.ics.forth.gr/~lourakis/levmar/. In order to compile using this library the HAVE_LEVMAR compiler flag must be set.

* The GNU Scientific Library (GSL) is also required and can be downloaded from https://www.gnu.org/software/gsl/. In order to compile using this library the HAVE_GSL compiler flag must be set.

## Intel Compiler
The binaries distributed on the FEBio website for OSX and Linux are compiled using the Intel compiler. This compiler can be downloaded from https://software.intel.com/en-us/c-compilers/. However, this compiler is not free. For instructions on how to compile FEBio using the GNU Compiler Collection (GCC), see the section below entitled "Linux and OSX Makefile Instructions".

## Windows Instructions
Included in the source are Visual Studio project files for Visual Studio 2008, 2010, 2013, and 2015. Most of the configuration required to build FEBio is set up in these project files. Depending on the location of the third party libraries on your system, it may be necessary to edit include and link paths in Visual Studio. Instructions on how to change these paths can be found on <a href="https://docs.microsoft.com/en-us/cpp/build/reference/vcpp-directories-property-page?view=vs-2019">Microsoft's Website</a>.

## OSX XCode Instructions
Included in the source is a project file for XCode for compilation on OSX. and 2015. Most of the configuration required to build FEBio is set up in the project file. Depending on the location of the third party libraries on your system, it may be necessary to edit include and link paths in XCode.

## Linux and OSX Makefile Instructions
Included under the build directory are makefile configurations for use with the GNU Make build system. All make configurations are defined in build/Makefile. Once invoked with a specifc configuration, this makefile calls several of the other makefiles in the build directory depending on the configuration that was specified.

Before any make configuration is called, a set of subdirectories must first be created in which the object files will be stored. This can be accomplished automatically by running the Mkdir.bash script located in the build directory. To do this open a terminal in the build directory and run the script with the name of the desired build configuration as an argument. For instance, in order to build the lnx64 configuration you would run:

```
./Mkdir.bash lnx64
```

Each configuration calls the `febio2.mk` and `Makelibs.mk` makefiles, which in turn include a configuration-specific makefile in which are defined include and link paths for the third pary libraries.

`lnx64d.mk` and `osxd.mk` are base configuration files for most of the make configurations on Linux and OSX respectively. For instance, when the lnx64 configuration is called, `lnx64.mk` simply includes the information in `lnx64d.mk` and removes the debug flag. This is done to make it easier to set up include and link paths. If the default include and link paths do not match the install locations of third party libraries on your machine, it  is only necessary to change these paths in `lnx64d.mk` in order for these changes to be made to all of the following configurations: lnx64, lnx64d, lnx64g, lnx64s, and gcc64.

A brief explanation of the available configurations follows:

lnx32:	Configuration for 32 bit Linux machines using the Intel compiler  
lnx64:	Configuration for 64 bit Linux machines using the Intel compiler  
gcc:	Configuration for 32 bit Linux machines using the GNU Compiler Collection (GCC)  
gcc64:	Configuration for 64 bit Linux machines using the GNU Compiler Collection (GCC)  
sky:	Configuration for 64 bit Linux machines using the GNU Compiler Collection (GCC), and the Skyline linear solver  

osx:	Configuration for 64 bit Macs using the Intel compiler  
osxg++:	Configuration for 64 bit Macs using the GNU Compiler Collection (GCC)  

The following suffixes may be appended to many of the previous configurations (e.g. lnx64s, osxd, etc):

d:	Configuration for a version of FEBio that allows FEBio to dump information about non-converged states to allow FE simulations to be debugged more easily. Note that this does not produce a debug executable.  
g:	Configuration used to produce a debug executable  
s:	Configuration used to build a sequential version of FEBio in which no multithreading is used.  

Once a build configuration has been decided on and the `Mkdir.bash` script has been run, simply run the make command followed by the configuration name as argument. For example, to build the lnx64 configuration you would run:

```
make lnx64
```

If your machine's processor has multiple cores, it is possible to decrease your build time by parallelizing your build with the -j flag followed by the number of cores you would like make to use. For example, to build the lnx64 configuration using 8 cores you would run:
```
make lnx64 -j8
```

Please note that this will only increase the speed of the compilation process and will in no way affect the final binary.

# Contributing <a name="Contributing"></a>

