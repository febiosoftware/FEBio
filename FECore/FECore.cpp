#include "stdafx.h"
#include "FECore.h"
#include "FECoordSysMap.h"
#include "FECoreKernel.h"
#include "BC.h"
#include "FEInitialCondition.h"
#include "FECorePlot.h"

#define FECORE_VERSION		0
#define FECORE_SUBVERSION	1

//-----------------------------------------------------------------------------
// Get the version info
void FECore::get_version(int& version, int& subversion)
{
	version = FECORE_VERSION;
	subversion = FECORE_SUBVERSION;
}

//-----------------------------------------------------------------------------
// get the version string
const char* FECore::get_version_string()
{
	static const char fecore_str[4] = {'0'+FECORE_VERSION, '.', '0'+FECORE_SUBVERSION };
	return fecore_str;
}

//-----------------------------------------------------------------------------
void FECore::InitModule()
{
REGISTER_FECORE_CLASS(FELocalMap         , FECOORDSYSMAP_ID, "local"      );
REGISTER_FECORE_CLASS(FESphericalMap     , FECOORDSYSMAP_ID, "spherical"  );
REGISTER_FECORE_CLASS(FECylindricalMap   , FECOORDSYSMAP_ID, "cylindrical");
REGISTER_FECORE_CLASS(FEVectorMap        , FECOORDSYSMAP_ID, "vector"     );
REGISTER_FECORE_CLASS(FESphericalAngleMap, FECOORDSYSMAP_ID, "angles"     );
REGISTER_FECORE_CLASS(FEPolarMap         , FECOORDSYSMAP_ID, "polar"      );

REGISTER_FECORE_CLASS(FEFixedBC      , FEBC_ID, "fix"      );
REGISTER_FECORE_CLASS(FEPrescribedDOF, FEBC_ID, "prescribe");
REGISTER_FECORE_CLASS(FENodalLoad    , FEBC_ID, "nodal load");

REGISTER_FECORE_CLASS(FEInitialBC     , FEIC_ID, "init_bc"       );
REGISTER_FECORE_CLASS(FEInitialBCVec3D, FEIC_ID, "init_bc_vec3d" );

REGISTER_FECORE_CLASS(FEPlotMaterialParameter, FEPLOTDATA_ID, "parameter");
}
