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
#include "FESolidModule.h"
#include <FECore/DOFS.h>
#include <FECore/FEModel.h>
#include "FEBioMech.h"

FESolidModule::FESolidModule() {}

void FESolidModule::InitModel(FEModel* fem)
{
	// Allocate degrees of freedom
	DOFS& dofs = fem->GetDOFS();
	int varD = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT), VAR_VEC3);
	dofs.SetDOFName(varD, 0, "x");
	dofs.SetDOFName(varD, 1, "y");
	dofs.SetDOFName(varD, 2, "z");
	int varQ = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_ROTATION), VAR_VEC3);
	dofs.SetDOFName(varQ, 0, "u");
	dofs.SetDOFName(varQ, 1, "v");
	dofs.SetDOFName(varQ, 2, "w");
	int varQR = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION), VAR_VEC3);
	dofs.SetDOFName(varQR, 0, "Ru");
	dofs.SetDOFName(varQR, 1, "Rv");
	dofs.SetDOFName(varQR, 2, "Rw");
	int varV = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::VELOCITY), VAR_VEC3);
	dofs.SetDOFName(varV, 0, "vx");
	dofs.SetDOFName(varV, 1, "vy");
	dofs.SetDOFName(varV, 2, "vz");
	int varSU = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT), VAR_VEC3);
	dofs.SetDOFName(varSU, 0, "sx");
	dofs.SetDOFName(varSU, 1, "sy");
	dofs.SetDOFName(varSU, 2, "sz");
	int varSV = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_VELOCITY), VAR_VEC3);
	dofs.SetDOFName(varSV, 0, "svx");
	dofs.SetDOFName(varSV, 1, "svy");
	dofs.SetDOFName(varSV, 2, "svz");
	int varSA = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_ACCELERATION), VAR_VEC3);
	dofs.SetDOFName(varSA, 0, "sax");
	dofs.SetDOFName(varSA, 1, "say");
	dofs.SetDOFName(varSA, 2, "saz");
	int varBW = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::BEAM_ANGULAR_VELOCITY), VAR_VEC3);
	dofs.SetDOFName(varBW, 0, "bwx");
	dofs.SetDOFName(varBW, 1, "bwy");
	dofs.SetDOFName(varBW, 2, "bwz");
	int varBA = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::BEAM_ANGULAR_ACCELERATION), VAR_VEC3);
	dofs.SetDOFName(varBA, 0, "bax");
	dofs.SetDOFName(varBA, 1, "bay");
	dofs.SetDOFName(varBA, 2, "baz");
}
