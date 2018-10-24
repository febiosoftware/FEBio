// Implementation of the Harwell-Boeing Matrix Print Diagnostic

#include "FEPrintHBMatrixDiagnostic.h"
#include "FEBioMech/FESolidSolver2.h"
#include <FECore/CompactMatrix.h>
#include <FECore/FEGlobalMatrix.h>

//-----------------------------------------------------------------------------
bool write_hb(CompactMatrix& K, const char* szfile)
{
	FILE* fp = fopen(szfile, "wb");
	if (fp == 0) return false;

	int	symmFlag = K.isSymmetric();
	int offset = K.Offset();
	int rowFlag = K.isRowBased();
	int neq = K.Rows();
	int nnz = K.NonZeroes();
	fwrite(&symmFlag, sizeof(symmFlag), 1, fp);
	fwrite(&offset, sizeof(offset), 1, fp);
	fwrite(&rowFlag, sizeof(rowFlag), 1, fp);
	fwrite(&neq, sizeof(neq), 1, fp);
	fwrite(&nnz, sizeof(nnz), 1, fp);
	fwrite(K.Pointers(), sizeof(int), neq + 1, fp);
	fwrite(K.Indices(), sizeof(int), nnz, fp);
	fwrite(K.Values(), sizeof(double), nnz, fp);

	fclose(fp);

	return true;
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
	FEModel& fem = GetFEModel();
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
	SparseMatrix* psm = solver->GetStiffnessMatrix().GetSparseMatrixPtr();
	CompactMatrix* pcm = dynamic_cast<CompactMatrix*>(psm);
	if (pcm == 0) return false;

	// print the matrix to file
	if (write_hb(*pcm, "hb_matrix.out") == false)
	{
		fprintf(stderr, "Failed writing sparse matrix.\n\n");
	}

	return true;
}
