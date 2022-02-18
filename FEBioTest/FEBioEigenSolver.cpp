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
#include "FEBioEigenSolver.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/EigenSolver.h>
#include <FECore/SparseMatrix.h>
#include <FEBioMech/FESolidSolver2.h>
#include <FEBioPlot/FEBioPlotFile.h>

FEBioEigenSolver::FEBioEigenSolver(FEModel* fem) : FECoreTask(fem)
{

}

bool FEBioEigenSolver::Init(const char* szfile)
{
	FEModel* fem = GetFEModel();
	if (fem == nullptr) return false;

	return fem->Init();
}

bool FEBioEigenSolver::Run()
{
	FEModel* fem = GetFEModel();
	if (fem == nullptr) return false;

	FEAnalysis* step = fem->GetStep(0);
	if (step->Activate() == false) return false;

	if (step->InitSolver() == false) return false;

	FESolidSolver2* solver = dynamic_cast<FESolidSolver2*>(fem->GetStep(0)->GetFESolver());
	if (solver == nullptr) return false;

	// evaluate to stiffness matrix
	if (solver->ReformStiffness() == false) return false;

	// get the stiffness matrix
	SparseMatrix* K = solver->GetStiffnessMatrix()->GetSparseMatrixPtr(); assert(K);

	// create the eigen solver
	EigenSolver* eigenSolver = fecore_new<EigenSolver>("feast", fem);
	if (eigenSolver == nullptr) return false;
	eigenSolver->GetParameter("m0")->value<int>() = K->Rows();
	eigenSolver->GetParameter("emin")->value<double>() = 0.0;
	eigenSolver->GetParameter("emax")->value<double>() = 1.0;

	// initialize eigen solver
	if (eigenSolver->Init() == false) return false;

	// get eigen values and eigen vectors
	vector<double> eigenValues;
	matrix eigenVectors;
	bool b = eigenSolver->EigenSolve(K, nullptr, eigenValues, eigenVectors);
	if (b == false) return false;

	FEMesh& mesh = fem->GetMesh();

	// write eigen values and eigen vectors
	FEBioPlotFile plt(fem);

	if (plt.AddVariable("displacement") == false) return false;

	if (plt.Open("eigen.xplt") == false) return false;

	plt.Write(0.0f);

	int n = -1;
	for (int i = 0; i < eigenValues.size(); ++i)
	{
		fem->GetTime().currentTime = eigenValues[i];

		for (int j = 0; j < mesh.Nodes(); ++j)
		{
			FENode& node = mesh.Node(j);
			node.m_rt = node.m_r0;
			n = node.m_ID[0]; if (n >= 0) node.m_rt.x += eigenVectors[i][n];
			n = node.m_ID[1]; if (n >= 0) node.m_rt.y += eigenVectors[i][n];
			n = node.m_ID[2]; if (n >= 0) node.m_rt.z += eigenVectors[i][n];
		}

		plt.Write(eigenValues[i]);
	}

	plt.Close();

	return true;
}
