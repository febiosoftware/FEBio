#include "stdafx.h"
#include "FECore/febio.h"

#include "FEBioPlot/FEPlotSurfaceData.h"
#include "FECore/FECoordSysMap.h"

#include "FEBioMech/FEBioMech.h"
#include "FEBioMix/FEBioMix.h"
#include "FEBioHeat/FEBioHeat.h"

void InitFEBioLibrary()
{
//-----------------------------------------------------------------------------
// import all modules
FEBioMech::InitModule();
FEBioMix::InitModule();
FEBioHeat::InitModule();

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FEPlotContactGap      , FEPlotData, "contact gap"     );
REGISTER_FEBIO_CLASS(FEPlotContactPressure , FEPlotData, "contact pressure");
REGISTER_FEBIO_CLASS(FEPlotContactTraction , FEPlotData, "contact traction");

//-----------------------------------------------------------------------------
// Classes derived from FECoordSysMap
REGISTER_FEBIO_CLASS(FELocalMap      , FECoordSysMap, "local"      );
REGISTER_FEBIO_CLASS(FESphericalMap  , FECoordSysMap, "spherical"  );
REGISTER_FEBIO_CLASS(FECylindricalMap, FECoordSysMap, "cylindrical");
REGISTER_FEBIO_CLASS(FEVectorMap     , FECoordSysMap, "vector"     );
}
