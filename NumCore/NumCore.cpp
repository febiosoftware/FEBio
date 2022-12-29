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
#include "NumCore.h"
#include "PardisoSolver.h"
#include "PardisoProjectSolver.h"
#include "RCICGSolver.h"
#include "FGMRESSolver.h"
#include "ILU0_Preconditioner.h"
#include "ILUT_Preconditioner.h"
#include "BIPNSolver.h"
#include "HypreGMRESsolver.h"
#include "Hypre_PCG_AMG.h"
#include "SchurSolver.h"
#include "IncompleteCholesky.h"
#include "BoomerAMGSolver.h"
#include "BlockSolver.h"
#include "BiCGStabSolver.h"
#include "StrategySolver.h"
#include <FECore/fecore_enum.h>
#include <FECore/FECoreFactory.h>
#include <FECore/FECoreKernel.h>
#include "FEASTEigenSolver.h"
#include "TestSolver.h"
#include "AccelerateSparseSolver.h"
#include "SuperLU_MT.h"
#include "MKLDSSolver.h"
#include "numcore_api.h"

//=============================================================================
// Call this to initialize the NumCore module
NUMCORE_API void NumCore::InitModule()
{
	// register linear solvers
	REGISTER_FECORE_CLASS(PardisoSolver  , "pardiso");
    REGISTER_FECORE_CLASS(PardisoProjectSolver, "pardiso-project");
	REGISTER_FECORE_CLASS(FGMRESSolver        , "fgmres"   );
	REGISTER_FECORE_CLASS(BoomerAMGSolver     , "boomeramg");
	REGISTER_FECORE_CLASS(RCICGSolver         , "cg"    );
	REGISTER_FECORE_CLASS(SchurSolver         , "schur"    );
	REGISTER_FECORE_CLASS(HypreGMRESsolver    , "hypre_gmres");
	REGISTER_FECORE_CLASS(Hypre_PCG_AMG       , "hypre_pcg_amg");
	REGISTER_FECORE_CLASS(BlockIterativeSolver, "block");
	REGISTER_FECORE_CLASS(BIPNSolver          , "bipn");
	REGISTER_FECORE_CLASS(BiCGStabSolver      , "bicgstab");
	REGISTER_FECORE_CLASS(StrategySolver      , "strategy");
	REGISTER_FECORE_CLASS(TestSolver          , "test");
    REGISTER_FECORE_CLASS(AccelerateSparseSolver, "accelerate");
    REGISTER_FECORE_CLASS(SuperLU_MT_Solver     , "superlu_mt");
    REGISTER_FECORE_CLASS(MKLDSSolver           , "mkl_dss");

	// register preconditioners
	REGISTER_FECORE_CLASS(ILU0_Preconditioner, "ilu0");
	REGISTER_FECORE_CLASS(ILUT_Preconditioner, "ilut");
	REGISTER_FECORE_CLASS(IncompleteCholesky , "ichol");

	// register eigen solvers
	REGISTER_FECORE_CLASS(FEASTEigenSolver, "feast");

	// set default linear solver
	// (Set this before the configuration is read in because
	//  the configuration can change the default linear solver.)
	FECoreKernel& fecore = FECoreKernel::GetInstance();
#ifdef PARDISO
	fecore.SetDefaultSolverType("pardiso");
#else
	fecore.SetDefaultSolverType("skyline");
#endif
}
