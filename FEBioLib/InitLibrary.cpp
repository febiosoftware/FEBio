#include "stdafx.h"
#include "FECore/FECore.h"
#include "FEBioMech/FEBioMech.h"
#include "FEBioMix/FEBioMix.h"
#include "FEBioHeat/FEBioHeat.h"
#include "FEMicroMaterial.h"

void InitFEBioLibrary()
{
//-----------------------------------------------------------------------------
// import all modules
FECore::InitModule();
FEBioMech::InitModule();
FEBioMix::InitModule();
FEBioHeat::InitModule();

//-----------------------------------------------------------------------------
// register the material with the framework
REGISTER_FEBIO_CLASS(FEMicroMaterial, FEMATERIAL_ID, "micro-material");
}
