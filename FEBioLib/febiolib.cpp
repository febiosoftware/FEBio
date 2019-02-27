#include "stdafx.h"
#include "FECore/FECore.h"
#include "NumCore/NumCore.h"
#include "FEBioMech/FEBioMech.h"
#ifndef MECH_ONLY
#include "FEBioMix/FEBioMix.h"
#include "FEBioOpt/FEBioOpt.h"
#include "FEBioFluid/FEBioFluid.h"
#include "FEBioFluid/FEBioFluidP.h"
#include <FEBioFluid/FEBioFSI.h>
#include <FEBioTest/FEBioTest.h>
#endif
#include "febio.h"
#include "plugin.h"
#include "FEBioStdSolver.h"

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
	REGISTER_FECORE_CLASS(FEBioRestart, "restart");

	FECore::InitModule();
	NumCore::InitModule();
	FEBioMech::InitModule();

#ifndef MECH_ONLY
	FEBioMix::InitModule();
	FEBioOpt::InitModule();
	FEBioFluid::InitModule();
    FEBioFluidP::InitModule();
	FEBioFSI::InitModule();
	FEBioTest::InitModule();
#endif
}

//-----------------------------------------------------------------------------
void FinishLibrary()
{
	FEBioPluginManager* pPM = FEBioPluginManager::GetInstance();
	pPM->DeleteThis();
}

} // namespace febio
