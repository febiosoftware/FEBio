#include "stdafx.h"
#include "FECore/FECore.h"
#include "FEBioMech/FEBioMech.h"
#include "FEBioMix/FEBioMix.h"
#include "FEBioHeat/FEBioHeat.h"
#include "FEBioOpt/FEBioOpt.h"
#include "FEBioFluid/FEBioFluid.h"

void InitFEBioLibrary()
{
//-----------------------------------------------------------------------------
// import all modules
FECore::InitModule();
FEBioMech::InitModule();
FEBioMix::InitModule();
FEBioHeat::InitModule();
FEBioOpt::InitModule();
FEBioFluid::InitModule();
}
