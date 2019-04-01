/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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

#include "stdafx.h"
#include "FEBioFSI.h"
#include <FECore/FECoreKernel.h>
#include "FEFluidFSISolver.h"
#include "FEFluidFSI.h"
#include "FEFluidFSIDomain3D.h"
#include "FEFluidFSITraction.h"
#include "FEFluidFSIDomainFactory.h"

void FEBioFSI::InitModule()
{
	FECoreKernel& febio = FECoreKernel::GetInstance();

	// register domain
	febio.RegisterDomain(new FEFluidFSIDomainFactory);

	// define the fsi module
	febio.CreateModule("fluid-FSI");
	febio.SetModuleDependency("fluid");

	REGISTER_FECORE_CLASS(FEFluidFSISolver, "fluid-FSI");

	REGISTER_FECORE_CLASS(FEFluidFSI, "fluid-FSI");

	REGISTER_FECORE_CLASS(FEFluidFSIDomain3D, "fluid-FSI-3D");

	REGISTER_FECORE_CLASS(FEFluidFSITraction, "fluid-FSI traction");

	febio.SetActiveModule(0);
}
