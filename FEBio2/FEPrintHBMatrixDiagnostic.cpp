// Implementation of the Harwell-Boeing Matrix Print Diagnostic

#include "stdafx.h"
#include "FEPrintHBMatrixDiagnostic.h"
#include "FEBioMech/FESolidSolver.h"
#include "NumCore/CompactMatrix.h"

//-----------------------------------------------------------------------------
void write_hb(CompactMatrix& m, FILE* fp)
{
	int neq = m.Size();
	int nnz = m.NonZeroes();

	fwrite(&neq, sizeof(neq), 1, fp);
	fwrite(&nnz, sizeof(nnz), 1, fp);
	fwrite(m.Pointers(), sizeof(int)   , neq+1, fp);
	fwrite(m.Indices (), sizeof(int)   , nnz, fp);
	fwrite(m.Values  (), sizeof(double), nnz, fp);
}

//-----------------------------------------------------------------------------
FEPrintHBMatrixDiagnostic::FEPrintHBMatrixDiagnostic(FEModel& fem) : FEDiagnostic(fem)
{
}

//-----------------------------------------------------------------------------
FEPrintHBMatrixDiagnostic::~FEPrintHBMatrixDiagnostic(void)
{
}

//-----------------------------------------------------------------------------
bool FEPrintHBMatrixDiagnostic::ParseSection(XMLTag &tag)
{
	FEModel& fem = m_fem;
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
			fprintf(stderr, szerr);

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
	FEAnalysis* pstep = m_fem.GetCurrentStep();

	// initialize the step
	pstep->Activate();

	// get and initialize the FE solver
	FESolidSolver& solver = static_cast<FESolidSolver&>(*pstep->m_psolver);
	solver.Init();

	// build the stiffness matrix
	// recalculate the shape of the stiffness matrix if necessary
	if (m_fem.SurfacePairInteractions()) solver.UpdateContact();

	// reshape the stiffness matrix
	if (!solver.CreateStiffness(true)) return false;

	// calculate the stiffness matrices
	FETimePoint tp; 
	tp.t = 0.0;
	tp.dt = 0.0;
	solver.StiffnessMatrix(tp);

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
