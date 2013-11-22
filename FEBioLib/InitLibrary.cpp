#include "stdafx.h"
#include "FECore/febio.h"

#include "FEBioPlot/FEPlotSurfaceData.h"

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
REGISTER_MATERIAL(FEMicroMaterial, "micro-material");

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FEPlotContactGap      , FEPlotData, "contact gap"     );
REGISTER_FEBIO_CLASS(FEPlotContactPressure , FEPlotData, "contact pressure");
REGISTER_FEBIO_CLASS(FEPlotContactTraction , FEPlotData, "contact traction");
}
