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
#include "FEBiphasicModule.h"
#include <FECore/DOFS.h>
#include <FECore/FEModel.h>
#include "FEBioMix.h"

//=============================================================================
FEBiphasicModule::FEBiphasicModule() {}

void FEBiphasicModule::InitModel(FEModel* fem)
{
	FESolidModule::InitModel(fem);

	// Allocate degrees of freedom
	DOFS& dofs = fem->GetDOFS();
	int varP = dofs.AddVariable("fluid pressure");
	dofs.SetDOFName(varP, 0, "p");
	int varQ = dofs.AddVariable("shell fluid pressure");
	dofs.SetDOFName(varQ, 0, "q");
}

//=============================================================================
FEBiphasicSoluteModule::FEBiphasicSoluteModule() {}

void FEBiphasicSoluteModule::InitModel(FEModel* fem)
{
	FEBiphasicModule::InitModel(fem);

	// Allocate degrees of freedom
	// (We start with zero concentration degrees of freedom)
	DOFS& dofs = fem->GetDOFS();
	int varC = dofs.AddVariable("concentration", VAR_ARRAY);
	int varD = dofs.AddVariable("shell concentration", VAR_ARRAY);
}

//=============================================================================
FEMultiphasicModule::FEMultiphasicModule() {}

void FEMultiphasicModule::InitModel(FEModel* fem)
{
	FEBiphasicSoluteModule::InitModel(fem);
}
