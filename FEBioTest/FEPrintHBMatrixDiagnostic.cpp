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
#include "FEPrintHBMatrixDiagnostic.h"
#include "FEBioMech/FESolidSolver2.h"
#include <FECore/CompactMatrix.h>
#include <FECore/FEGlobalMatrix.h>
#include <NumCore/MatrixTools.h>

//-----------------------------------------------------------------------------
FEPrintHBMatrixDiagnostic::FEPrintHBMatrixDiagnostic(FEModel* fem) : FEDiagnostic(fem)
{
}

//-----------------------------------------------------------------------------
FEPrintHBMatrixDiagnostic::~FEPrintHBMatrixDiagnostic(void)
{
}

//-----------------------------------------------------------------------------
bool FEPrintHBMatrixDiagnostic::ParseSection(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	if (tag == "input")
	{
		// get the input file name
		const char* szfile = tag.szvalue();

		// try to read the file
		FEBioImport im;
		if (im.Load(fem, szfile) == false)
		{
			char szerr[256];
			im.GetErrorMessage(szerr);
			fprintf(stderr, "%s", szerr);

			return false;
		}

		return true;
	}

	return  false;
}

//-----------------------------------------------------------------------------
bool FEPrintHBMatrixDiagnostic::Run()
{
	// Get the current step
	FEModel& fem = *GetFEModel();
	FEAnalysis* pstep = fem.GetCurrentStep();

	// initialize the step
	pstep->Activate();

	// get and initialize the FE solver
	FENewtonSolver* solver = dynamic_cast<FENewtonSolver*>(pstep->GetFESolver());
	if (solver == 0) return false;

	// do initialization
	FETimeInfo& tp = fem.GetTime();
	tp.currentTime = pstep->m_dt0;
	tp.timeIncrement = pstep->m_dt0;
	solver->InitStep(pstep->m_dt0);
	solver->PrepStep();

	// reshape the stiffness matrix
	if (!solver->ReformStiffness()) return false;

	// get the matrix
	SparseMatrix* psm = solver->GetStiffnessMatrix()->GetSparseMatrixPtr();
	CompactMatrix* pcm = dynamic_cast<CompactMatrix*>(psm);
	if (pcm == 0) return false;

	// print the matrix to file
	if (NumCore::write_hb(*pcm, "hb_matrix.out") == false)
	{
		fprintf(stderr, "Failed writing sparse matrix.\n\n");
	}

	return true;
}
