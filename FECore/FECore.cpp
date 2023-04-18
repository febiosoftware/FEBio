/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FECore.h"
#include "FECoreKernel.h"
#include "FEPrescribedDOF.h"
#include "FENodalLoad.h"
#include "FEFixedBC.h"
#include "FELinearConstraint.h"
#include "FEInitialCondition.h"
#include "FECorePlot.h"
#include "FESurfaceToSurfaceMap.h"
#include "FEParabolicMap.h"
#include "FEDataMathGenerator.h"
#include "FEPointFunction.h"
#include "FELoadCurve.h"
#include "FEMathController.h"
#include "FEMathIntervalController.h"
#include "FEPIDController.h"
#include "Preconditioner.h"
#include "FEMat3dValuator.h"
#include "FEMat3dSphericalAngleMap.h"
#include "FEAnalysis.h"
#include "BFGSSolver.h"
#include "FEBroydenStrategy.h"
#include "JFNKStrategy.h"
#include "FENodeSet.h"
#include "FEFacetSet.h"
#include "FEElementSet.h"
#include "FEConstValueVec3.h"
#include "NodeDataRecord.h"
#include "FaceDataRecord.h"
#include "ElementDataRecord.h"
#include "NLConstraintDataRecord.h"
#include "FEAugLagLinearConstraint.h"
#include "SurfaceDataRecord.h"
#include "FELogEnclosedVolume.h"
#include "FELogElementVolume.h"
#include "FELogDomainVolume.h"
#include "FELogSolutionNorm.h"
#include "FELinearConstraint.h"
#include "LUSolver.h"
#include "FETimeStepController.h"
#include "FEModifiedNewtonStrategy.h"
#include "FEFullNewtonStrategy.h"
#include "SkylineSolver.h"

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
//REGISTER_FECORE_CLASS(FEAnalysis, "analysis");

// time controller
REGISTER_FECORE_CLASS(FETimeStepController, "default");

// boundary conditions
REGISTER_FECORE_CLASS(FEFixedDOF     , "fix"      , 0x300);	// obsolete in 4.0
REGISTER_FECORE_CLASS(FEPrescribedDOF, "prescribe", 0x300);	// obsolete in 4.0
REGISTER_FECORE_CLASS(FELinearConstraint, "linear constraint");
REGISTER_FECORE_CLASS(FELinearConstraintDOF, "child_dof");

// nodal loads
REGISTER_FECORE_CLASS(FENodalDOFLoad, "nodal_load");

// initial conditions
REGISTER_FECORE_CLASS(FEInitialDOF     , "init_dof"     , 0x300);	// obsolete in 4.0

// (augmented lagrangian) linear constraints
REGISTER_FECORE_CLASS(FELinearConstraintSet, "linear constraint");
REGISTER_FECORE_CLASS(FEAugLagLinearConstraint, "linear_constraint");
REGISTER_FECORE_CLASS(FEAugLagLinearConstraintDOF, "node");

// plot field
REGISTER_FECORE_CLASS(FEPlotParameter, "parameter");

// load curves
REGISTER_FECORE_CLASS(FEPointFunction , "point");
REGISTER_FECORE_CLASS(FEConstFunction, "const");
REGISTER_FECORE_CLASS(FELinearFunction, "linear ramp");
REGISTER_FECORE_CLASS(FEStepFunction  , "step");
REGISTER_FECORE_CLASS(FEMathFunction  , "math");

// data generators
REGISTER_FECORE_CLASS(FEDataMathGenerator  , "math");
REGISTER_FECORE_CLASS(FESurfaceToSurfaceMap, "surface-to-surface map");
REGISTER_FECORE_CLASS(FEParabolicMap       , "parabolic map");

// scalar valuators
REGISTER_FECORE_CLASS(FEConstValue , "const");
REGISTER_FECORE_CLASS(FEMathValue  , "math" );
REGISTER_FECORE_CLASS(FEMappedValue, "map"  );

//  vector generators
REGISTER_FECORE_CLASS(FELocalVectorGenerator          , "local");
REGISTER_FECORE_CLASS(FEConstValueVec3                , "vector");
REGISTER_FECORE_CLASS(FEMathValueVec3                 , "math");
REGISTER_FECORE_CLASS(FESphericalVectorGenerator      , "spherical");
REGISTER_FECORE_CLASS(FECylindricalVectorGenerator    , "cylindrical");
REGISTER_FECORE_CLASS(FESphericalAnglesVectorGenerator, "angles");
REGISTER_FECORE_CLASS(FEMappedValueVec3               , "map");
REGISTER_FECORE_CLASS(FEUserVectorGenerator           , "user");

// mat3d generators
REGISTER_FECORE_CLASS(FEConstValueMat3d       , "const"      );
REGISTER_FECORE_CLASS(FEMat3dLocalElementMap  , "local"      );
REGISTER_FECORE_CLASS(FEMat3dSphericalMap     , "spherical"  );
REGISTER_FECORE_CLASS(FEMat3dCylindricalMap   , "cylindrical");
REGISTER_FECORE_CLASS(FEMat3dVectorMap        , "vector"     );
REGISTER_FECORE_CLASS(FEMat3dSphericalAngleMap, "angles"     );
REGISTER_FECORE_CLASS(FEMat3dPolarMap         , "polar"      );
REGISTER_FECORE_CLASS(FEMappedValueMat3d      , "map"        );

// mat3ds generators
REGISTER_FECORE_CLASS(FEConstValueMat3ds , "const");
REGISTER_FECORE_CLASS(FEMappedValueMat3ds, "map");

// load controllers
REGISTER_FECORE_CLASS(FELoadCurve             , "loadcurve");
REGISTER_FECORE_CLASS(FEMathController        , "math");
REGISTER_FECORE_CLASS(FEMathIntervalController, "math-interval");
REGISTER_FECORE_CLASS(FEPIDController         , "PID");

// Newton strategies
REGISTER_FECORE_CLASS(BFGSSolver       , "BFGS");
REGISTER_FECORE_CLASS(FEBroydenStrategy, "Broyden");
REGISTER_FECORE_CLASS(JFNKStrategy     , "JFNK");
REGISTER_FECORE_CLASS(FEModifiedNewtonStrategy, "modified Newton");
REGISTER_FECORE_CLASS(FEFullNewtonStrategy    , "full Newton");

// preconditioners
REGISTER_FECORE_CLASS(DiagonalPreconditioner, "diagonal");

REGISTER_FECORE_CLASS(FESurface, "surface");

// data records
REGISTER_FECORE_CLASS(NodeDataRecord, "node_data");
REGISTER_FECORE_CLASS(FaceDataRecord, "face_data");
REGISTER_FECORE_CLASS(ElementDataRecord, "element_data");
REGISTER_FECORE_CLASS(NLConstraintDataRecord, "rigid_connector_data");

// log classes
REGISTER_FECORE_CLASS(FELogEnclosedVolume, "volume");
REGISTER_FECORE_CLASS(FELogElementVolume, "V");
REGISTER_FECORE_CLASS(FELogDomainVolume, "volume");
REGISTER_FECORE_CLASS(FELogAvgDomainData, "avg");
REGISTER_FECORE_CLASS(FELogPctDomainData, "pct");
REGISTER_FECORE_CLASS(FELogIntegralDomainData, "integrate");
REGISTER_FECORE_CLASS(FELogSolutionNorm, "solution_norm");
REGISTER_FECORE_CLASS(FELogFaceArea    , "facet area");

// linear solvers
REGISTER_FECORE_CLASS(LUSolver, "LU");
REGISTER_FECORE_CLASS(SkylineSolver, "skyline");

}
