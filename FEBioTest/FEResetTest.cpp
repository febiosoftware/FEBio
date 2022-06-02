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
#include "FEResetTest.h"
#include <FEBioLib/FEBioModel.h>
#include <FECore/log.h>
#include <FEBioLib/Logfile.h>
#include <FECore/FELogSolutionNorm.h>
#include <iostream>
#include <iomanip>
using namespace std;

//-----------------------------------------------------------------------------
FEResetTest::FEResetTest(FEModel*pfem) : FECoreTask(pfem)
{
}

//-----------------------------------------------------------------------------
// initialize the diagnostic
bool FEResetTest::Init(const char* sz)
{
	FEBioModel& fem = dynamic_cast<FEBioModel&>(*GetFEModel());

	Logfile& log = fem.GetLogFile();
	log.SetMode(Logfile::MODE::LOG_FILE);

	// do the FE initialization
	return fem.Init();
}

//-----------------------------------------------------------------------------
// run the diagnostic
bool FEResetTest::Run()
{
	FEBioModel* fem = dynamic_cast<FEBioModel*>(GetFEModel());

	FELogSolutionNorm sn(fem);

	// try to run the model
	cerr << "Running model for the first time.\n";
	if (fem->Solve() == false)
	{
		feLogEx(fem, "Failed to run model.");
		return false;
	}

	// collect results
	ModelStats stats1 = fem->GetModelStats();
	double norm1 = sn.value();
	cerr << "time steps    = " << stats1.ntimeSteps    << endl;
	cerr << "total iters   = " << stats1.ntotalIters   << endl;
	cerr << "total reforms = " << stats1.ntotalReforms << endl;
	cerr << "total rhs     = " << stats1.ntotalRHS     << endl;
	cerr << "solution norm = " << std::setprecision(15) << norm1 << endl;

	// reset the model
	std::cerr << "Resetting model.\n";
	if (fem->Reset() == false)
	{
		feLogEx(fem, "Failed to reset model.");
		return false;
	}

	// try to run it again
	std::cerr << "Running model for the second time.\n";
	if (fem->Solve() == false)
	{
		feLogEx(fem, "Failed to run model second time.");
		return false;
	}

	// get model stats
	ModelStats stats2 = fem->GetModelStats();
	double norm2 = sn.value();
	cerr << "time steps    = " << stats2.ntimeSteps    << endl;
	cerr << "total iters   = " << stats2.ntotalIters   << endl;
	cerr << "total reforms = " << stats2.ntotalReforms << endl;
	cerr << "total rhs     = " << stats2.ntotalRHS     << endl;
	cerr << "solution norm = " << norm2 << endl;

	bool success = true;
	if (stats1.ntimeSteps    != stats2.ntimeSteps   ) success = false;
	if (stats1.ntotalIters   != stats2.ntotalIters  ) success = false;
	if (stats1.ntotalReforms != stats2.ntotalReforms) success = false;
	if (stats1.ntotalRHS     != stats2.ntotalRHS    ) success = false;
	if (norm1 != norm2) success = false;
	cerr << " --> Reset test " << (success ? "PASSED" : "FAILED") << endl;

	return success;
}
