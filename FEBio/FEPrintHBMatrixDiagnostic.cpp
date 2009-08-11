#include "stdafx.h"
#include "FEPrintHBMatrixDiagnostic.h"

FEPrintHBMatrixDiagnostic::FEPrintHBMatrixDiagnostic(FEM& fem) : FEDiagnostic(fem)
{
}

FEPrintHBMatrixDiagnostic::~FEPrintHBMatrixDiagnostic(void)
{
}

bool FEPrintHBMatrixDiagnostic::ParseSection(XMLTag &tag)
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

	return  false;
}

bool FEPrintHBMatrixDiagnostic::Run()
{
	// get and initialize the first step
	m_fem.m_Step[0].Init();

	// get and initialize the solver
	FESolver& solver = *m_fem.m_pStep->m_psolver;
	solver.Init();

	// build the stiffness matrix
	// recalculate the shape of the stiffness matrix if necessary
	if (m_fem.m_bcontact) m_fem.UpdateContact();

	// reshape the stiffness matrix
	if (!solver.CreateStiffness(true)) return false;

	// calculate the stiffness matrices
	solver.StiffnessMatrix();

	// print the matrix
	SparseMatrix* psm = solver.GetStiffnessMatrix()->GetSparseMatrixPtr();
	CompactMatrix* pcm = dynamic_cast<CompactMatrix*>(psm);
	if (pcm == 0) return false;
	if (!pcm->print_hb()) return false;

	return true;
}
