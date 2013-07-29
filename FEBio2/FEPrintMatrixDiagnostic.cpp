#include "stdafx.h"
#include "FEPrintMatrixDiagnostic.h"
#include "FEBioMech/FESolidSolver.h"

FEPrintMatrixDiagnostic::FEPrintMatrixDiagnostic(FEModel& fem) : FEDiagnostic(fem)
{
	m_szout[0] = 0;
	m_rng[0] = m_rng[1] = 0;
	m_rng[2] = m_rng[3] = -1;
}

FEPrintMatrixDiagnostic::~FEPrintMatrixDiagnostic(void)
{
}

bool FEPrintMatrixDiagnostic::ParseSection(XMLTag &tag)
{
	if (tag == "input")
	{
		// get the input file name
		const char* szfile = tag.szvalue();

		// try to read the file
		FEFEBioImport im;
		if (im.Load(m_fem, szfile) == false)
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

bool FEPrintMatrixDiagnostic::Run()
{
	// get and initialize the first step
	FEAnalysisStep* pstep = dynamic_cast<FEAnalysisStep*>(m_fem.GetStep(0));
	pstep->Init();

	// get and initialize the solver
	FESolidSolver& solver = dynamic_cast<FESolidSolver&>(*pstep->m_psolver);
	solver.Init();

	// build the stiffness matrix
	// recalculate the shape of the stiffness matrix if necessary
	if (m_fem.SurfacePairInteractions() > 0) solver.UpdateContact();

	// reshape the stiffness matrix
	if (!solver.CreateStiffness(true)) return false;

	// calculate the stiffness matrices
	FETimePoint tp = {0.0, 0.0};
	solver.StiffnessMatrix(tp);

	// print the matrix
	FILE* fp = fopen(m_szout, "wt");
	if (fp == 0) { fprintf(stderr, "Failed creating output file."); return false; }
	int* n = m_rng;
	print(*solver.GetStiffnessMatrix()->GetSparseMatrixPtr(), fp, n[0], n[1], n[2], n[3]);
	fclose(fp);

	return true;
}
