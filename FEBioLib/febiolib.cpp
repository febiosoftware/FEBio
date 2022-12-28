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
#include "FECore/FECore.h"
#include "NumCore/NumCore.h"
#include "FEAMR/FEAMR.h"
#include "FEBioMech/FEBioMechModule.h"
#ifndef MECH_ONLY
#include "FEBioMix/FEBioMix.h"
#include "FEBioOpt/FEBioOpt.h"
#include "FEBioFluid/FEBioFluid.h"
#include <FEBioFluid/FEBioFSI.h>
#include <FEBioFluid/FEBioMultiphasicFSI.h>
#include <FEBioFluid/FEBioFluidSolutes.h>
#include <FEBioFluid/FEBioThermoFluid.h>
#include <FEBioFluid/FEBioPolarFluid.h>
#include <FEBioTest/FEBioTest.h>
#include <FEBioRVE/FEBioRVE.h>
#include <FEImgLib/FEImgLib.h>
#endif
#include "febio.h"
#include "plugin.h"
#include "FEBioStdSolver.h"
#include "FEBioRestart.h"

namespace febio {

//-----------------------------------------------------------------------------
FECoreKernel* GetFECoreKernel()
{
	return &FECoreKernel::GetInstance();
}

//-----------------------------------------------------------------------------
// import all modules
void InitLibrary()
{
	REGISTER_FECORE_CLASS(FEBioStdSolver, "solve");
	REGISTER_FECORE_CLASS(FEBioRestart  , "restart");
	REGISTER_FECORE_CLASS(FEBioRCISolver, "rci_solve");
	REGISTER_FECORE_CLASS(FEBioTestSuiteTask, "test");

	FECore::InitModule();
	FEAMR::InitModule();
	NumCore::InitModule();
	FEBioMech::InitModule();

#ifndef MECH_ONLY
	FEBioMix::InitModule();
	FEBioOpt::InitModule();
	FEBioFluid::InitModule();
	FEBioFSI::InitModule();
    FEBioMultiphasicFSI::InitModule();
    FEBioFluidSolutes::InitModule();
    FEBioThermoFluid::InitModule();
    FEBioPolarFluid::InitModule();
	FEBioTest::InitModule();
	FEBioRVE::InitModule();
	FEImgLib::InitModule();
#endif
}

//-----------------------------------------------------------------------------
void FinishLibrary()
{
	FEBioPluginManager* pPM = FEBioPluginManager::GetInstance();
	pPM->DeleteThis();
}

} // namespace febio
