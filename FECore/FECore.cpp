#include "stdafx.h"
#include "FECore.h"
#include "FECoordSysMap.h"
#include "FECoreKernel.h"
#include "FEPrescribedDOF.h"
#include "FENodalLoad.h"
#include "FEFixedBC.h"
#include "FEInitialCondition.h"
#include "FECorePlot.h"
#include "FEDataLoadCurve.h"

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
// coordinate system map
REGISTER_FECORE_CLASS(FELocalMap         , FECOORDSYSMAP_ID, "local"      );
REGISTER_FECORE_CLASS(FESphericalMap     , FECOORDSYSMAP_ID, "spherical"  );
REGISTER_FECORE_CLASS(FECylindricalMap   , FECOORDSYSMAP_ID, "cylindrical");
REGISTER_FECORE_CLASS(FEVectorMap        , FECOORDSYSMAP_ID, "vector"     );
REGISTER_FECORE_CLASS(FESphericalAngleMap, FECOORDSYSMAP_ID, "angles"     );
REGISTER_FECORE_CLASS(FEPolarMap         , FECOORDSYSMAP_ID, "polar"      );

// boundary conditions
REGISTER_FECORE_CLASS(FEFixedBC      , FEBC_ID, "fix"      );
REGISTER_FECORE_CLASS(FEPrescribedDOF, FEBC_ID, "prescribe");
REGISTER_FECORE_CLASS(FENodalLoad    , FEBC_ID, "nodal load");

// initial conditions
REGISTER_FECORE_CLASS(FEInitialBC     , FEIC_ID, "init_bc"       );
REGISTER_FECORE_CLASS(FEInitialBCVec3D, FEIC_ID, "init_bc_vec3d" );

// plot field
REGISTER_FECORE_CLASS(FEPlotParameter, FEPLOTDATA_ID, "parameter");

// load curves
REGISTER_FECORE_CLASS(FEDataLoadCurve, FELOADCURVE_ID, "loadcurve");
REGISTER_FECORE_CLASS(FELinearRamp   , FELOADCURVE_ID, "linear ramp");
}
