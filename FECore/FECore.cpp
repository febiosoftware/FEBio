/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "BFGSSolver.h"
#include "FEBroydenStrategy.h"
#include "JFNKStrategy.h"
#include "FENodeSet.h"
#include "FEFacetSet.h"
#include "FEElementSet.h"
#include "FETetgenRefine.h"

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

// nodal loads
REGISTER_FECORE_CLASS(FENodalDOFLoad, "nodal_load");

// initial conditions
REGISTER_FECORE_CLASS(FEInitialDOF     , "init_dof"     );

// plot field
REGISTER_FECORE_CLASS(FEPlotParameter, "parameter");

// load curves
REGISTER_FECORE_CLASS(FEPointFunction , "point");
REGISTER_FECORE_CLASS(FELinearFunction, "linear ramp");

// data generators
REGISTER_FECORE_CLASS(FEDataMathGenerator  , "math");
REGISTER_FECORE_CLASS(FESurfaceToSurfaceMap, "surface-to-surface map");

// scalar valuators
REGISTER_FECORE_CLASS(FEConstValue , "const");
REGISTER_FECORE_CLASS(FEMathValue  , "math" );
REGISTER_FECORE_CLASS(FEMappedValue, "map"  );

//  vector generators
REGISTER_FECORE_CLASS(FELocalVectorGenerator      , "local");
REGISTER_FECORE_CLASS(FEConstValueVec3            , "vector");
REGISTER_FECORE_CLASS(FEMathValueVec3             , "math");
REGISTER_FECORE_CLASS(FESphericalVectorGenerator  , "spherical");
REGISTER_FECORE_CLASS(FECylindricalVectorGenerator, "cylindrical");
REGISTER_FECORE_CLASS(FEMappedValueVec3           , "map");

// mat3d generators
REGISTER_FECORE_CLASS(FEConstValueMat3d       , "const"      );
REGISTER_FECORE_CLASS(FEMat3dLocalElementMap  , "local"      );
REGISTER_FECORE_CLASS(FEMat3dSphericalMap     , "spherical"  );
REGISTER_FECORE_CLASS(FEMat3dCylindricalMap   , "cylindrical");
REGISTER_FECORE_CLASS(FEMat3dVectorMap        , "vector"     );
REGISTER_FECORE_CLASS(FEMat3dSphericalAngleMap, "angles"     );
REGISTER_FECORE_CLASS(FEMat3dPolarMap         , "polar"      );
REGISTER_FECORE_CLASS(FEMappedValueMat3d      , "map"        );

// load controllers
REGISTER_FECORE_CLASS(FELoadCurve     , "loadcurve");
REGISTER_FECORE_CLASS(FEMathController, "math");
REGISTER_FECORE_CLASS(FEPIDController , "PID");

// Newton strategies
REGISTER_FECORE_CLASS(BFGSSolver       , "BFGS");
REGISTER_FECORE_CLASS(FEBroydenStrategy, "Broyden");
REGISTER_FECORE_CLASS(JFNKStrategy     , "JFNK");

// preconditioners
REGISTER_FECORE_CLASS(DiagonalPreconditioner, "diagonal");

// Mesh item lists
REGISTER_FECORE_CLASS(FENodeSet   , "node_set");
REGISTER_FECORE_CLASS(FEFacetSet  , "surface" );
REGISTER_FECORE_CLASS(FEElementSet, "elem_set");

// mesh adaptors
REGISTER_FECORE_CLASS(FEErosionAdaptor, "erosion");
REGISTER_FECORE_CLASS(FEHexRefine     , "hex_refine");
REGISTER_FECORE_CLASS(FEHexRefine2D   , "hex_refine2d");
REGISTER_FECORE_CLASS(FETetRefine     , "tet_refine");
REGISTER_FECORE_CLASS(FETetgenRefine  , "tetgen_refine");

REGISTER_FECORE_CLASS(FEMaxVolumeCriterion, "max_volume");
REGISTER_FECORE_CLASS(FEMaxVariableCriterion, "max_variable");
}
