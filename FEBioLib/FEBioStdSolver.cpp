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
#include "FEBioStdSolver.h"
#include <FEBioLib/FEBioModel.h>
#include <FECore/log.h>
#include <FEBioXML/FERestartImport.h>
#include <FECore/DumpFile.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEModelDataRecord.h>

//-----------------------------------------------------------------------------
FEBioStdSolver::FEBioStdSolver(FEModel* pfem) : FECoreTask(pfem) {}

//-----------------------------------------------------------------------------
// This simply calls the FEModel::Init
bool FEBioStdSolver::Init(const char* szfile)
{
	return (GetFEModel() ? GetFEModel()->Init() : false);
}

//-----------------------------------------------------------------------------
// This simply calls the FEM::Solve function which will solve the FE problem.
bool FEBioStdSolver::Run()
{
	// Solve the problem and return error code
	return (GetFEModel() ? GetFEModel()->Solve() : false);
}

//=============================================================================

FEBioRCISolver::FEBioRCISolver(FEModel* fem) : FECoreTask(fem) {}

//! initialization
bool FEBioRCISolver::Init(const char* szfile)
{
	return (GetFEModel() ? GetFEModel()->Init() : false);
}

//! Run the FE model
bool FEBioRCISolver::Run()
{
	// get the model
	FEModel* fem = GetFEModel();
	if (fem == nullptr) return false;

	// initialize RCI solver
	if (fem->RCI_Init() == false) return false;

	// loop until solved
	while (fem->IsSolved() == false)
	{
		// try to advance the solution
		if (fem->RCI_Advance() == false)
		{
			// if we were unable to advance the solution, we do a rewind and try again
			if (fem->RCI_Rewind() == false)
			{
				// couldn't rewind, so we're done
				break;
			}
		}
	}

	// finalize the solver
	if (fem->RCI_Finish() == false) return false;

	return true;
}

//==========================================================================
FEBioTestSuiteTask::FEBioTestSuiteTask(FEModel* fem) : FECoreTask(fem) {}

//! initialization
bool FEBioTestSuiteTask::Init(const char* szfile)
{
	FEModel* fem = GetFEModel(); assert(fem);
	if (fem == nullptr) return false;

	// See if the model defines any data records
	DataStore& data = fem->GetDataStore();
	if (data.Size() == 0)
	{
		FEModelDataRecord* rec = new FEModelDataRecord(fem);
		rec->SetData("solution_norm");
		rec->SetName("solution_norm");
		data.AddRecord(rec);
	}

	return (GetFEModel() ? GetFEModel()->Init() : false);
}

//! Run the FE model
bool FEBioTestSuiteTask::Run()
{
	return (GetFEModel() ? GetFEModel()->Solve() : false);
}
