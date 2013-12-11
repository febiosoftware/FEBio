#include "stdafx.h"
#include "FECore/FECore.h"
#include "FEBioMech/FEBioMech.h"
#include "FEBioMix/FEBioMix.h"
#include "FEBioHeat/FEBioHeat.h"


void InitFEBioLibrary()
{
//-----------------------------------------------------------------------------
// import all modules
FECore::InitModule();
FEBioMech::InitModule();
FEBioMix::InitModule();
FEBioHeat::InitModule();
}
