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
#include "FEMemoryDiagnostic.h"
#include "FECore/log.h"

FEMemoryDiagnostic::FEMemoryDiagnostic(FEModel* fem) : FEDiagnostic(fem)
{
	m_szfile[0] = 0;
	m_iters = 1;

	FEAnalysis* pstep = new FEAnalysis(fem);
	fem->AddStep(pstep);
	fem->SetCurrentStep(pstep);
}

FEMemoryDiagnostic::~FEMemoryDiagnostic(void)
{

}

bool FEMemoryDiagnostic::Init()
{
	FEModel& fem = *GetFEModel();

	// try to open the file
	FEBioImport im;
	if (im.Load(fem, m_szfile) == false)
	{
		return false;
	}

	// make sure the iters is a positive number
	if (m_iters <= 0) return false;

	// turn off all output
	fem.GetCurrentStep()->SetPlotLevel(FE_PLOT_NEVER);

	return true;
}

bool FEMemoryDiagnostic::Run()
{
	// run the problem
	FEModel& fem = *GetFEModel();
	for (int i=0; i<m_iters; ++i)
	{
		fprintf(stderr, "%d/%d: ...", i+1, m_iters);
		fem.Reset();
		bool b = fem.Solve();
//		system("ps -C febio.test -o rss,vsize h");
		fprintf(stderr, "%s\n", (b?"NT" : "ET"));
	}

	return true;
}

bool FEMemoryDiagnostic::ParseSection(XMLTag &tag)
{
	if (tag == "file") tag.value(m_szfile);
	else if (tag == "iters") tag.value(m_iters);
	else return false;
	return true;
}
