/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#include "FEFluidModule.h"
#include <FECore/FEModel.h>
#include "FEBioFluid.h"
#include "FEBioFSI.h"
#include "FEBioMultiphasicFSI.h"
#include "FEBioThermoFluid.h"
#include "FEBioFluidSolutes.h"
#include "FEBioPolarFluid.h"

FEFluidModule::FEFluidModule() {}
void FEFluidModule::InitModel(FEModel* fem)
{
    // Allocate degrees of freedom
    DOFS& dofs = fem->GetDOFS();
    int varD = dofs.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::DISPLACEMENT), VAR_VEC3);
    dofs.SetDOFName(varD, 0, "x");
    dofs.SetDOFName(varD, 1, "y");
    dofs.SetDOFName(varD, 2, "z");

    int nW = dofs.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(nW, 0, "wx");
    dofs.SetDOFName(nW, 1, "wy");
    dofs.SetDOFName(nW, 2, "wz");

    int nE = dofs.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION), VAR_SCALAR);
    dofs.SetDOFName(nE, 0, "ef");

    int nAW = dofs.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(nAW, 0, "awx");
    dofs.SetDOFName(nAW, 1, "awy");
    dofs.SetDOFName(nAW, 2, "awz");

    int nAE = dofs.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION_TDERIV), VAR_SCALAR);
    dofs.SetDOFName(nAE, 0, "aef");
}

//=============================================================================
FEFluidFSIModule::FEFluidFSIModule() {}
void FEFluidFSIModule::InitModel(FEModel* fem)
{
    // Allocate degrees of freedom
    DOFS& dofs = fem->GetDOFS();

    int varD = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::DISPLACEMENT), VAR_VEC3);
    dofs.SetDOFName(varD, 0, "x");
    dofs.SetDOFName(varD, 1, "y");
    dofs.SetDOFName(varD, 2, "z");

    int varV = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::VELOCITY), VAR_VEC3);
    dofs.SetDOFName(varV, 0, "vx");
    dofs.SetDOFName(varV, 1, "vy");
    dofs.SetDOFName(varV, 2, "vz");

    int varQ = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_ROTATION), VAR_VEC3);
    dofs.SetDOFName(varQ, 0, "u");
    dofs.SetDOFName(varQ, 1, "v");
    dofs.SetDOFName(varQ, 2, "w");

    int varSD = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_DISPLACEMENT), VAR_VEC3);
    dofs.SetDOFName(varSD, 0, "sx");
    dofs.SetDOFName(varSD, 1, "sy");
    dofs.SetDOFName(varSD, 2, "sz");

    int varQV = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(varQV, 0, "svx");
    dofs.SetDOFName(varQV, 1, "svy");
    dofs.SetDOFName(varQV, 2, "svz");

    int varQA = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(varQA, 0, "sax");
    dofs.SetDOFName(varQA, 1, "say");
    dofs.SetDOFName(varQA, 2, "saz");

    int varQR = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RIGID_ROTATION), VAR_VEC3);
    dofs.SetDOFName(varQR, 0, "Ru");
    dofs.SetDOFName(varQR, 1, "Rv");
    dofs.SetDOFName(varQR, 2, "Rw");

    int nW = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RELATIVE_FLUID_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(nW, 0, "wx");
    dofs.SetDOFName(nW, 1, "wy");
    dofs.SetDOFName(nW, 2, "wz");

    int nAW = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RELATIVE_FLUID_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(nAW, 0, "awx");
    dofs.SetDOFName(nAW, 1, "awy");
    dofs.SetDOFName(nAW, 2, "awz");

    int nVF = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::FLUID_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(nVF, 0, "vfx");
    dofs.SetDOFName(nVF, 1, "vfy");
    dofs.SetDOFName(nVF, 2, "vfz");

    int nAF = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::FLUID_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(nAF, 0, "afx");
    dofs.SetDOFName(nAF, 1, "afy");
    dofs.SetDOFName(nAF, 2, "afz");

    int nE = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::FLUID_DILATATION), VAR_SCALAR);
    dofs.SetDOFName(nE, 0, "ef");

    int nAE = dofs.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::FLUID_DILATATION_TDERIV), VAR_SCALAR);
    dofs.SetDOFName(nAE, 0, "aef");
}

//=============================================================================
FEMultiphasicFSIModule::FEMultiphasicFSIModule() { SetStatus(EXPERIMENTAL); }
void FEMultiphasicFSIModule::InitModel(FEModel* fem)
{
    // Allocate degrees of freedom
    DOFS& dofs = fem->GetDOFS();

    int varD = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::DISPLACEMENT), VAR_VEC3);
    dofs.SetDOFName(varD, 0, "x");
    dofs.SetDOFName(varD, 1, "y");
    dofs.SetDOFName(varD, 2, "z");

    int varV = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::VELOCITY), VAR_VEC3);
    dofs.SetDOFName(varV, 0, "vx");
    dofs.SetDOFName(varV, 1, "vy");
    dofs.SetDOFName(varV, 2, "vz");

    int varQ = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::SHELL_ROTATION), VAR_VEC3);
    dofs.SetDOFName(varQ, 0, "u");
    dofs.SetDOFName(varQ, 1, "v");
    dofs.SetDOFName(varQ, 2, "w");

    int varSD = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::SHELL_DISPLACEMENT), VAR_VEC3);
    dofs.SetDOFName(varSD, 0, "sx");
    dofs.SetDOFName(varSD, 1, "sy");
    dofs.SetDOFName(varSD, 2, "sz");

    int varQV = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::SHELL_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(varQV, 0, "svx");
    dofs.SetDOFName(varQV, 1, "svy");
    dofs.SetDOFName(varQV, 2, "svz");

    int varQA = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::SHELL_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(varQA, 0, "sax");
    dofs.SetDOFName(varQA, 1, "say");
    dofs.SetDOFName(varQA, 2, "saz");

    int varQR = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::RIGID_ROTATION), VAR_VEC3);
    dofs.SetDOFName(varQR, 0, "Ru");
    dofs.SetDOFName(varQR, 1, "Rv");
    dofs.SetDOFName(varQR, 2, "Rw");

    int nW = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::RELATIVE_FLUID_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(nW, 0, "wx");
    dofs.SetDOFName(nW, 1, "wy");
    dofs.SetDOFName(nW, 2, "wz");

    int nAW = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::RELATIVE_FLUID_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(nAW, 0, "awx");
    dofs.SetDOFName(nAW, 1, "awy");
    dofs.SetDOFName(nAW, 2, "awz");

    int nVF = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(nVF, 0, "vfx");
    dofs.SetDOFName(nVF, 1, "vfy");
    dofs.SetDOFName(nVF, 2, "vfz");

    int nAF = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(nAF, 0, "afx");
    dofs.SetDOFName(nAF, 1, "afy");
    dofs.SetDOFName(nAF, 2, "afz");

    int nE = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_DILATATION), VAR_SCALAR);
    dofs.SetDOFName(nE, 0, "ef");

    int nAE = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_DILATATION_TDERIV), VAR_SCALAR);
    dofs.SetDOFName(nAE, 0, "aef");

    int varC = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_CONCENTRATION), VAR_ARRAY);
    int varAC = dofs.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_CONCENTRATION_TDERIV), VAR_ARRAY);
}

//=============================================================================
FEThermoFluidModule::FEThermoFluidModule() { SetStatus(EXPERIMENTAL); }
void FEThermoFluidModule::InitModel(FEModel* fem)
{
    // Allocate degrees of freedom
    DOFS& dofs = fem->GetDOFS();
    int varD = dofs.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::DISPLACEMENT), VAR_VEC3);
    dofs.SetDOFName(varD, 0, "x");
    dofs.SetDOFName(varD, 1, "y");
    dofs.SetDOFName(varD, 2, "z");

    int nW = dofs.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::RELATIVE_FLUID_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(nW, 0, "wx");
    dofs.SetDOFName(nW, 1, "wy");
    dofs.SetDOFName(nW, 2, "wz");

    int nE = dofs.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::FLUID_DILATATION), VAR_SCALAR);
    dofs.SetDOFName(nE, 0, "ef");

    int nT = dofs.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE), VAR_SCALAR);
    dofs.SetDOFName(nT, 0, "T");

    int nAW = dofs.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::RELATIVE_FLUID_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(nAW, 0, "awx");
    dofs.SetDOFName(nAW, 1, "awy");
    dofs.SetDOFName(nAW, 2, "awz");

    int nAE = dofs.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::FLUID_DILATATION_TDERIV), VAR_SCALAR);
    dofs.SetDOFName(nAE, 0, "aef");

    int nAT = dofs.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE_TDERIV), VAR_SCALAR);
    dofs.SetDOFName(nAT, 0, "aT");
}

//=============================================================================
FEPolarFluidModule::FEPolarFluidModule() { SetStatus(EXPERIMENTAL); }
void FEPolarFluidModule::InitModel(FEModel* fem)
{
    // Allocate degrees of freedom
    DOFS& dofs = fem->GetDOFS();
    
    int varD = dofs.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::DISPLACEMENT), VAR_VEC3);
    dofs.SetDOFName(varD, 0, "x");
    dofs.SetDOFName(varD, 1, "y");
    dofs.SetDOFName(varD, 2, "z");
    
    int nW = dofs.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::RELATIVE_FLUID_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(nW, 0, "wx");
    dofs.SetDOFName(nW, 1, "wy");
    dofs.SetDOFName(nW, 2, "wz");
    
    int nAW = dofs.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::RELATIVE_FLUID_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(nAW, 0, "awx");
    dofs.SetDOFName(nAW, 1, "awy");
    dofs.SetDOFName(nAW, 2, "awz");
    
    int nG = dofs.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::FLUID_ANGULAR_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(nG, 0, "gx");
    dofs.SetDOFName(nG, 1, "gy");
    dofs.SetDOFName(nG, 2, "gz");
    
    int nAG = dofs.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::FLUID_ANGULAR_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(nAG, 0, "agx");
    dofs.SetDOFName(nAG, 1, "agy");
    dofs.SetDOFName(nAG, 2, "agz");
    
    int nE = dofs.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::FLUID_DILATATION), VAR_SCALAR);
    dofs.SetDOFName(nE, 0, "ef");
    
    int nAE = dofs.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::FLUID_DILATATION_TDERIV), VAR_SCALAR);
    dofs.SetDOFName(nAE, 0, "aef");
}

//=============================================================================
FEFluidSolutesModule::FEFluidSolutesModule() { SetStatus(EXPERIMENTAL); }
//FEFluidSolutesModule::FEFluidSolutesModule() { SetStatus(RELEASED); }
void FEFluidSolutesModule::InitModel(FEModel* fem)
{
    // Allocate degrees of freedom
    DOFS& dofs = fem->GetDOFS();
    int varD = dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::DISPLACEMENT), VAR_VEC3);
    dofs.SetDOFName(varD, 0, "x");
    dofs.SetDOFName(varD, 1, "y");
    dofs.SetDOFName(varD, 2, "z");

    int nW = dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::RELATIVE_FLUID_VELOCITY), VAR_VEC3);
    dofs.SetDOFName(nW, 0, "wx");
    dofs.SetDOFName(nW, 1, "wy");
    dofs.SetDOFName(nW, 2, "wz");

    int nE = dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_DILATATION), VAR_SCALAR);
    dofs.SetDOFName(nE, 0, "ef");

    int nAW = dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::RELATIVE_FLUID_ACCELERATION), VAR_VEC3);
    dofs.SetDOFName(nAW, 0, "awx");
    dofs.SetDOFName(nAW, 1, "awy");
    dofs.SetDOFName(nAW, 2, "awz");

    int nAE = dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_DILATATION_TDERIV), VAR_SCALAR);
    dofs.SetDOFName(nAE, 0, "aef");

    int varC = dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), VAR_ARRAY);
    int varAC = dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION_TDERIV), VAR_ARRAY);
}
