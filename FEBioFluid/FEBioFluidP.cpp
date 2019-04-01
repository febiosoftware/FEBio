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
#include "FEBioFluidP.h"
#include <FECore/FECoreKernel.h>
#include "FEFluidSolver.h"
#include "FEFluidPDomain3D.h"
#include "FEFluidPDomainFactory.h"
#include "FEFluidPResistanceBC.h"

void FEBioFluidP::InitModule()
{
    FECoreKernel& febio = FECoreKernel::GetInstance();
    
    // register domain
    febio.RegisterDomain(new FEFluidPDomainFactory);
    
    // define the fluidP module
    febio.CreateModule("fluidP");
    febio.SetModuleDependency("fluid");

    REGISTER_FECORE_CLASS(FEFluidPDomain3D, "fluidP-3D");
    
    REGISTER_FECORE_CLASS(FEFluidPResistanceBC, "fluidP resistance");
    
    febio.SetActiveModule(0);
}
