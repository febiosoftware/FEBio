// Implementation of the Harwell-Boeing Matrix Print Diagnostic

#include "stdafx.h"
#include "FEPrintHBMatrixDiagnostic.h"
#include "FESolidSolver.h"
#include "FECore/CompactMatrix.h"
using namespace FECore;

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
	m_fem.m_Step[0]->Init();

	// get and initialize the solver
	FESolidSolver& solver = dynamic_cast<FESolidSolver&>(*m_fem.m_pStep->m_psolver);
	solver.Init();

	// build the stiffness matrix
	// recalculate the shape of the stiffness matrix if necessary
	if (m_fem.ContactInterfaces()) m_fem.UpdateContact();

	// reshape the stiffness matrix
	if (!solver.CreateStiffness(true)) return false;

	// calculate the stiffness matrices
	solver.StiffnessMatrix();

	// get the matrix
	SparseMatrix* psm = solver.GetStiffnessMatrix()->GetSparseMatrixPtr();
	CompactSymmMatrix* pcm = dynamic_cast<CompactSymmMatrix*>(psm);
	if (pcm == 0) return false;

	// print the matrix to file
	FILE* fout = fopen("hb_matrix.out", "wb");
	if (fout == 0) { fprintf(stderr, "Failed creating output file."); return false; }
	write_hb(*pcm, fout);
	fclose(fout);

	return true;
}
