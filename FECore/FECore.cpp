#include "stdafx.h"
#include "FECore.h"
#include "FECoordSysMap.h"
#include "FECoreKernel.h"
#include "BC.h"

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

REGISTER_FECORE_CLASS(FEFixedBC     , FEBC_ID, "fix"      );
REGISTER_FECORE_CLASS(FEPrescribedBC, FEBC_ID, "prescribe");
REGISTER_FECORE_CLASS(FENodalLoad   , FEBC_ID, "nodal load");
}
