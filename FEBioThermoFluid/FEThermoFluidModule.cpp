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
#include "FEThermoFluidModule.h"
#include "FEBioThermoFluid.h"
#include <FECore/FEModel.h>

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
