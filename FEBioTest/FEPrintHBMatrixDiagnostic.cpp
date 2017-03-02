// Implementation of the Harwell-Boeing Matrix Print Diagnostic

#include "FEPrintHBMatrixDiagnostic.h"
#include "FEBioMech/FESolidSolver2.h"
#include "NumCore/CompactMatrix.h"
#include "FECore/FEGlobalMatrix.h"

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
	FEAnalysis* pstep = new FEAnalysis(&fem);
    fem.AddStep(pstep);
    fem.SetCurrentStep(pstep);
}

//-----------------------------------------------------------------------------
FEPrintHBMatrixDiagnostic::~FEPrintHBMatrixDiagnostic(void)
{
}

//-----------------------------------------------------------------------------
bool FEPrintHBMatrixDiagnostic::ParseSection(XMLTag &tag)
{
	FEModel& fem = GetFEModel();
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
	FEModel& fem = GetFEModel();
	FEAnalysis* pstep = fem.GetCurrentStep();

	// initialize the step
	pstep->Init();
	pstep->Activate();

	// get and initialize the FE solver
	FESolidSolver2& solver = static_cast<FESolidSolver2&>(*pstep->GetFESolver());
	solver.Init();

	// build the stiffness matrix
	// recalculate the shape of the stiffness matrix if necessary
	if (fem.SurfacePairInteractions()) solver.UpdateContact();

	// reshape the stiffness matrix
	if (!solver.CreateStiffness(true)) return false;

	// calculate the stiffness matrices
	FETimeInfo tp; 
	solver.StiffnessMatrix(tp);

	// get the matrix
	SparseMatrix* psm = solver.GetStiffnessMatrix().GetSparseMatrixPtr();
	CompactSymmMatrix* pcm = dynamic_cast<CompactSymmMatrix*>(psm);
	if (pcm == 0) return false;

	// print the matrix to file
	FILE* fout = fopen("hb_matrix.out", "wb");
	if (fout == 0) { fprintf(stderr, "Failed creating output file."); return false; }
	write_hb(*pcm, fout);
	fclose(fout);

	return true;
}
