#include "stdafx.h"
#include "FECore.h"
#include "FECoreKernel.h"
#include "FEPrescribedDOF.h"
#include "FENodalLoad.h"
#include "FEFixedBC.h"
#include "FEInitialCondition.h"
#include "FECorePlot.h"
#include "FESurfaceToSurfaceMap.h"
#include "FEDataMathGenerator.h"
#include "FEPointFunction.h"
#include "FELoadCurve.h"
#include "FEMathController.h"
#include "FEPIDController.h"
#include "Preconditioner.h"
#include "FEMat3dValuator.h"
#include "FEAnalysis.h"
#include "FEErosionAdaptor.h"
#include "FEHexRefine.h"
#include "FEHexRefine2D.h"
#include "FETetRefine.h"

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
	// initialize the element librar
	FEElementLibrary::Initialize();

// analysis class
REGISTER_FECORE_CLASS(FEAnalysis, "analysis");

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

// scalar valuators
REGISTER_FECORE_CLASS(FEConstValue, "const");
REGISTER_FECORE_CLASS(FEMathValue , "math");

//  vector generators
REGISTER_FECORE_CLASS(FELocalVectorGenerator      , "local");
REGISTER_FECORE_CLASS(FEConstValueVec3            , "vector");
REGISTER_FECORE_CLASS(FEMathValueVec3             , "math");
REGISTER_FECORE_CLASS(FESphericalVectorGenerator  , "spherical");
REGISTER_FECORE_CLASS(FECylindricalVectorGenerator, "cylindrical");
REGISTER_FECORE_CLASS(FEMappedValueVec3           , "user");

// mat3d generators
REGISTER_FECORE_CLASS(FEConstValueMat3d       , "const"      );
REGISTER_FECORE_CLASS(FEMat3dLocalElementMap  , "local"      );
REGISTER_FECORE_CLASS(FEMat3dSphericalMap     , "spherical"  );
REGISTER_FECORE_CLASS(FEMat3dCylindricalMap   , "cylindrical");
REGISTER_FECORE_CLASS(FEMat3dVectorMap        , "vector"     );
REGISTER_FECORE_CLASS(FEMat3dSphericalAngleMap, "angles"     );
REGISTER_FECORE_CLASS(FEMat3dPolarMap         , "polar"      );
REGISTER_FECORE_CLASS(FEMappedValueMat3d      , "user"       );

// load controllers
REGISTER_FECORE_CLASS(FELoadCurve     , "loadcurve");
REGISTER_FECORE_CLASS(FEMathController, "math");
REGISTER_FECORE_CLASS(FEPIDController , "PID");

// preconditioners
REGISTER_FECORE_CLASS(DiagonalPreconditioner, "diagonal");

// mesh adaptors
REGISTER_FECORE_CLASS(FEErosionAdaptor, "erosion");
REGISTER_FECORE_CLASS(FEHexRefine     , "hex_refine");
REGISTER_FECORE_CLASS(FEHexRefine2D   , "hex_refine2d");
REGISTER_FECORE_CLASS(FETetRefine     , "tet_refine");

REGISTER_FECORE_CLASS(FEMaxVolumeCriterion, "max_volume");
}
