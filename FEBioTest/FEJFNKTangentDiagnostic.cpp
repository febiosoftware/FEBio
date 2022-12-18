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
#include "FEJFNKTangentDiagnostic.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FENewtonSolver.h>
#include <FECore/JFNKMatrix.h>
#include <FECore/CompactUnSymmMatrix.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FEGlobalMatrix.h>
#include <XML/XMLReader.h>

FEJFNKTangentDiagnostic::FEJFNKTangentDiagnostic(FEModel* fem) : FECoreTask(fem)
{
	m_max_err = 0.0;
	m_sum_err = 0.0;
	m_max_err_K = 0.0;
	m_max_err_A = 0.0;
	m_max_err_i = -1;
	m_max_err_j = -1;
	m_evals = 0;
	m_max_K = 0.0;

	m_jfnk_eps = 1e-7;
	m_small_val = 1e-9;
	m_nstep = -1;

	m_nprint = 0;
}

bool FEJFNKTangentDiagnostic::Init(const char* szfile)
{
	XMLReader xml;
	if (xml.Open(szfile) == false)
	{
		fprintf(stderr, "\nERROR: Failed to open %s\n\n", szfile);
		return false;
	}

	XMLTag tag;
	if (xml.FindTag("jfnk_diagnostic_spec", tag) == false)
	{
		fprintf(stderr, "\nERROR: Failed to read %s\n\n", szfile);
		return false;
	}

	++tag;
	do
	{
		if      (tag == "jfnk_eps" ) tag.value(m_jfnk_eps);
		else if (tag == "small_value") tag.value(m_small_val);
		else if (tag == "time_step") tag.value(m_nstep);
		else if (tag == "print_level") tag.value(m_nprint);
		else
		{
			fprintf(stderr, "ERROR: Failed to read %s\n\n", szfile);
			return false;
		}
		++tag;
	}
	while (!tag.isend());

	printf("JFNK Diagnostic parameters:\n");
	printf("\tJFNK eps    = %lg\n", m_jfnk_eps);
	printf("\tSmall value = %lg\n", m_small_val);
	printf("\ttime step   = %d\n", m_nstep);
	printf("\n");

	xml.Close();

	return GetFEModel()->Init();
}

bool cb_diagnose(FEModel* fem, unsigned int nwhen, void* pd)
{
	return ((FEJFNKTangentDiagnostic*)pd)->Diagnose();
}

bool FEJFNKTangentDiagnostic::Run()
{
	FEModel& fem = *GetFEModel();
	fem.AddCallback(cb_diagnose, CB_MAJOR_ITERS, this);

//	fem.GetLogFile().SetMode(Logfile::LOG_FILE);
	printf("Running model ...");
	bool ret = fem.Solve();
	printf("done!\n");

	printf("\nmax K    : %lg\n", m_max_K);

	printf("\nmax error: %lg\n", m_max_err);
	printf("\tK_exact[%d, %d] = %lg\n", m_max_err_i, m_max_err_j, m_max_err_K);
	printf("\tK_jfnk [%d, %d] = %lg\n", m_max_err_i, m_max_err_j, m_max_err_A);

	double avg_err = m_sum_err / m_evals;
	printf("\navg error: %lg\n\n", avg_err);

	return true;
}

bool FEJFNKTangentDiagnostic::Diagnose()
{
	FEModel& fem = *GetFEModel();

	// only analyze the requested step
	int timeSteps = fem.GetCurrentStep()->m_ntimesteps;
	if ((m_nstep != -1) && (timeSteps != m_nstep)) return true;

	// get the solver
	FENewtonSolver* fesolver = dynamic_cast<FENewtonSolver*>(fem.GetCurrentStep()->GetFESolver());
	if (fesolver == nullptr) return false;

	// Force a reform of the stiffness matrix first
	fesolver->ReformStiffness();

	// get the actual stiffness matrix
	FEGlobalMatrix* G = fesolver->GetStiffnessMatrix();
	if (G == nullptr) return false;

	SparseMatrix* K = G->GetSparseMatrixPtr();
	if (K == nullptr) return false;

	// setup the JFNK matrix
	JFNKMatrix A(fesolver);
	A.SetEpsilon(m_jfnk_eps);

	// get number of equations
	int neq = K->Rows();
	vector<double> ej(neq, 0.0), aj(neq, 0.0), kj(neq, 0.0);
	vector<bool> cj(neq, false);

	if (m_nprint)
	{
		FILE* fp = fopen("out_exact.txt", "wt");
		// print the transpose of the matrix 
		for (int j = 0; j < neq; ++j)
		{
			for (int i = 0; i < neq; ++i)
			{
				double kij = K->get(i, j);
				fprintf(fp, "%13.7lg", kij);
				if (i != neq - 1) fprintf(fp, ",");
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
	}

	FILE* out = NULL;
	if (m_nprint)
	{
		out = fopen("out_jfnk.txt", "wt");
	}

	// re-evaluate residual
	vector<double> R(neq, 0.0);
	fesolver->Residual(R);
	A.SetReferenceResidual(R);

	// loop over all columns
	for (int j = 0; j < neq; ++j)
	{
		// create e_i vector
		ej[j] = 1.0;

		// multiply with JFNK matrix
		A.mult_vector(&ej[0], &aj[0]);

		// get column of K
		for (int i = 0; i < neq; ++i)
		{
			cj[i] = K->check(i, j);
			double kij = K->get(i, j);
			kj[i] = kij;
			if (fabs(kij) > m_max_K) m_max_K = fabs(kij);
		}

		// compare difference
		for (int i = 0; i < neq; ++i)
		{
			double aij = aj[i];
			double kij = kj[i];

			if (m_nprint)
			{
				fprintf(out, "%13.7lg", aij);
				if (i != neq - 1) fprintf(out, ",");
			}

			double eij = 0;
			if (fabs(kij) <= m_small_val)
			{
				eij = fabs(kij - aij);
			}
			else 
			{
				double D = 0.5*(fabs(kij) + fabs(aij));
				if (D != 0.0)
				{
					eij = fabs(kij - aij) / D;
				}
			}

			if (eij != 0.0)
			{
				m_sum_err += eij;
				m_evals++;

				if (eij > m_max_err)
				{
					m_max_err = eij;
					m_max_err_K = kij;
					m_max_err_A = aij;
					m_max_err_i = i;
					m_max_err_j = j;
				}
			}
		}

		if (m_nprint)
		{
			fprintf(out, "\n");
		}
		
		// reset e_j vector
		ej[j] = 0.0;
	}

	if (m_nprint) fclose(out);

	return true;
}
