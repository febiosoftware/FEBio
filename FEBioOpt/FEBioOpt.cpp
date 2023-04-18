/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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



#include "stdafx.h"
#include "FEBioOpt.h"
#include "FEOptimize.h"
#include "FEParameterSweep.h"
#include "FEBioParamRun.h"
#include <FECore/FECoreKernel.h>
#include "FELMOptimizeMethod.h"
#include "FEConstrainedLMOptimizeMethod.h"
#include "FEPowellOptimizeMethod.h"
#include "FEScanOptimizeMethod.h"

//-----------------------------------------------------------------------------
//! Initialization of the FEBioOpt module. This function registers all the classes
//! in this module with the FEBio framework.
void FEBioOpt::InitModule()
{
	// Task classes
	REGISTER_FECORE_CLASS(FEOptimize      , "optimize");
	REGISTER_FECORE_CLASS(FEParameterSweep, "parameter_sweep");
	REGISTER_FECORE_CLASS(FEBioParamRun   , "param_run");

	// optimization methods
	REGISTER_FECORE_CLASS(FELMOptimizeMethod, "levmar");
#ifdef HAVE_LEVMAR
	REGISTER_FECORE_CLASS(FEConstrainedLMOptimizeMethod, "constrained levmar");
#endif
	REGISTER_FECORE_CLASS(FEPowellOptimizeMethod, "powell");
	REGISTER_FECORE_CLASS(FEScanOptimizeMethod, "scan");
}
