#include "FEPrintMatrixDiagnostic.h"
#include "FEBioMech/FESolidSolver2.h"
#include "FECore/FEGlobalMatrix.h"
#include <stdio.h>

//-----------------------------------------------------------------------------
//! Print a block from a sparse matrix to file
void print(SparseMatrix& m, FILE* fp, int i0, int j0, int i1, int j1)
{
	int ndim = m.Size();
	if ((i1 < 0) || (i1 >= ndim)) i1 = ndim-1;
	if ((j1 < 0) || (j1 >= ndim)) j1 = ndim-1;

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
FEPrintMatrixDiagnostic::FEPrintMatrixDiagnostic(FEModel& fem) : FEDiagnostic(fem)
{
	m_szout[0] = 0;
	m_rng[0] = m_rng[1] = 0;
	m_rng[2] = m_rng[3] = -1;

	FEAnalysis* pstep = new FEAnalysis(&fem);
    fem.AddStep(pstep);
    fem.m_nStep = 0;
    fem.SetCurrentStep(pstep);
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
		FEModel& fem = GetFEModel();
		if (im.Load(fem, szfile) == false)
		{
			char szerr[256];
			im.GetErrorMessage(szerr);
			fprintf(stderr, szerr);

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
	FEModel& fem = GetFEModel();
	FEAnalysis* pstep = fem.GetStep(0);
	pstep->Init();
	pstep->Activate();

	// get and initialize the solver
	FESolidSolver2& solver = static_cast<FESolidSolver2&>(*pstep->GetFESolver());
	solver.Init();

	// build the stiffness matrix
	// recalculate the shape of the stiffness matrix if necessary
	if (fem.SurfacePairInteractions() > 0) solver.UpdateContact();

	// reshape the stiffness matrix
	if (!solver.CreateStiffness(true)) return false;

	// calculate the stiffness matrices
	FETimePoint tp; 
	tp.t = 0.0;
	tp.dt = 0.0;
	solver.StiffnessMatrix(tp);

	// print the matrix
	FILE* fp = fopen(m_szout, "wt");
	if (fp == 0) { fprintf(stderr, "Failed creating output file."); return false; }
	int* n = m_rng;
	print(*solver.GetStiffnessMatrix().GetSparseMatrixPtr(), fp, n[0], n[1], n[2], n[3]);
	fclose(fp);

	return true;
}
