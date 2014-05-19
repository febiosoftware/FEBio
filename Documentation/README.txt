************************************************ FEBio README *************************************

FEBio is a nonlinear finite element solver that is specifically designed for biomechanical applications. It offers modeling scenarios, constitutive models and boundary conditions that are relevant to many research areas in biomechanics. All features can be used together seamlessly, giving the user a powerful tool for solving 3D problems in computational biomechanics. The software is open-source, and pre-compiled executables for Windows, OS-X and Linux platforms are available.

FEBio can be downloaded from http://mrl.sci.utah.edu/software/febio.  Please inform us of publications that use FEBio in research.  Information can be found on the Publications tab.  Support forums are at http://mrlforums.sci.utah.edu/forums/.


FEBio Makefile Guide
********************

FEBio comes with the default linear solver Skyline.  Instructions are included in build/Makefile to add the Pardiso linear solver.  The source code download also contains Visual Studio and XCode configurations.

febio.xml
*********

FEBio now requires a configuration file febio.xml.  This needs to be in the same directory as the executable or can be specified on the command line such as: febio.lnx -i inputfile.feb -cnf someconfig.xml.  The febio.xml file in the bin directory specifies Pardiso as the linear solver.  This can be edited to specify another linear solver.  For example, replace <linear_solver type="pardiso"></linear_solver> with <linear_solver type="skyline"></linear_solver> to use the skyline solver.