#include "stdafx.h"
#include "FECore.h"
#include "FECoordSysMap.h"
#include "FECoreKernel.h"
#include "FEPrescribedDOF.h"
#include "FENodalLoad.h"
#include "FEFixedBC.h"
#include "FEInitialCondition.h"
#include "FECorePlot.h"
#include "FESurfaceToSurfaceMap.h"
#include "FEDataMathGenerator.h"
#include "FEPointFunction.h"
#include "FEVectorGenerator.h"
#include "FELoadCurve.h"
#include "FEMathController.h"
#include "FEPIDController.h"
#include "Preconditioner.h"

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
REGISTER_FECORE_CLASS(FELocalMap         , "local"      );
REGISTER_FECORE_CLASS(FESphericalMap     , "spherical"  );
REGISTER_FECORE_CLASS(FECylindricalMap   , "cylindrical");
REGISTER_FECORE_CLASS(FEVectorMap        , "vector"     );
REGISTER_FECORE_CLASS(FESphericalAngleMap, "angles"     );
REGISTER_FECORE_CLASS(FEPolarMap         , "polar"      );

// boundary conditions
REGISTER_FECORE_CLASS(FEFixedBC      , "fix"      );
REGISTER_FECORE_CLASS(FEPrescribedDOF, "prescribe");
REGISTER_FECORE_CLASS(FENodalLoad    , "nodal load");

// initial conditions
REGISTER_FECORE_CLASS(FEInitialBC     , "init_bc"       );
REGISTER_FECORE_CLASS(FEInitialBCVec3D, "init_bc_vec3d" );

// plot field
REGISTER_FECORE_CLASS(FEPlotParameter, "parameter");

// load curves
REGISTER_FECORE_CLASS(FEPointFunction , "point");
REGISTER_FECORE_CLASS(FELinearFunction, "linear ramp");

// data generators
REGISTER_FECORE_CLASS(FEDataMathGenerator  , "math");
REGISTER_FECORE_CLASS(FESurfaceToSurfaceMap, "surface-to-surface map");

//  vector generators
REGISTER_FECORE_CLASS(FELocalVectorGenerator      , "local");
REGISTER_FECORE_CLASS(FEConstVectorGenerator      , "vector");
REGISTER_FECORE_CLASS(FESphericalVectorGenerator  , "spherical");
REGISTER_FECORE_CLASS(FECylindricalVectorGenerator, "cylindrical");
REGISTER_FECORE_CLASS(FEUserVectorGenerator       , "user");

// load controllers
REGISTER_FECORE_CLASS(FELoadCurve     , "loadcurve");
REGISTER_FECORE_CLASS(FEMathController, "math");
REGISTER_FECORE_CLASS(FEPIDController , "PID");

// preconditioners
REGISTER_FECORE_CLASS(DiagonalPreconditioner, "diagonal");
}
