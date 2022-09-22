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
#include "FEPrintMatrixDiagnostic.h"
#include "FEBioMech/FESolidSolver2.h"
#include "FECore/FEGlobalMatrix.h"
#include <stdio.h>

//-----------------------------------------------------------------------------
//! Print a block from a sparse matrix to file
void print(SparseMatrix& m, FILE* fp, int i0, int j0, int i1, int j1)
{
	int nr = m.Rows();
	int nc = m.Columns();
	if ((i1 < 0) || (i1 >= nr)) i1 = nr-1;
	if ((j1 < 0) || (j1 >= nc)) j1 = nc-1;

	for (int i=i0; i<=i1; ++i)
	{
		for (int j=j0; j<=j1; ++j)
		{
			fprintf(fp, "%10.3g", m.get(i,j));
		}
		fprintf(fp, "\n");
	}
}

//-----------------------------------------------------------------------------
FEPrintMatrixDiagnostic::FEPrintMatrixDiagnostic(FEModel* fem) : FEDiagnostic(fem)
{
	m_szout[0] = 0;
	m_rng[0] = m_rng[1] = 0;
	m_rng[2] = m_rng[3] = -1;

	FEAnalysis* pstep = new FEAnalysis(fem);
	fem->AddStep(pstep);
	fem->SetCurrentStep(pstep);
}

//-----------------------------------------------------------------------------
FEPrintMatrixDiagnostic::~FEPrintMatrixDiagnostic(void)
{
}

//-----------------------------------------------------------------------------
bool FEPrintMatrixDiagnostic::ParseSection(XMLTag &tag)
{
	if (tag == "input")
	{
		// get the input file name
		const char* szfile = tag.szvalue();

		// try to read the file
		FEBioImport im;
		FEModel& fem = *GetFEModel();
		if (im.Load(fem, szfile) == false)
		{
			char szerr[256];
			im.GetErrorMessage(szerr);
			fprintf(stderr, "%s", szerr);

			return false;
		}

		return true;
	}
	else if (tag == "output")
	{
		strcpy(m_szout, tag.szvalue());
		return true;
	}
	else if (tag == "range")
	{
		tag.value(m_rng, 4);
		return true;
	}

	return  false;
}

//-----------------------------------------------------------------------------
bool FEPrintMatrixDiagnostic::Run()
{
	// get and initialize the first step
	FEModel& fem = *GetFEModel();
	FEAnalysis* pstep = fem.GetStep(0);
	pstep->Init();
	pstep->Activate();

	// get and initialize the solver
	FESolidSolver2& solver = static_cast<FESolidSolver2&>(*pstep->GetFESolver());
	solver.Init();

	// build the stiffness matrix
	// recalculate the shape of the stiffness matrix if necessary
	fem.Update();

	// reshape the stiffness matrix
	if (!solver.CreateStiffness(true)) return false;

	// calculate the stiffness matrices
	solver.StiffnessMatrix();

	// print the matrix
	FILE* fp = fopen(m_szout, "wt");
	if (fp == 0) { fprintf(stderr, "Failed creating output file."); return false; }
	int* n = m_rng;
	print(*solver.GetStiffnessMatrix()->GetSparseMatrixPtr(), fp, n[0], n[1], n[2], n[3]);
	fclose(fp);

	return true;
}
