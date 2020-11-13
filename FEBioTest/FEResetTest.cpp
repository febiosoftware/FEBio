/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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

//-----------------------------------------------------------------------------
FEResetTest::FEResetTest(FEModel*pfem) : FECoreTask(pfem)
{
}

//-----------------------------------------------------------------------------
// initialize the diagnostic
bool FEResetTest::Init(const char* sz)
{
	FEBioModel& fem = dynamic_cast<FEBioModel&>(*GetFEModel());

	// do the FE initialization
	return fem.Init();
}

//-----------------------------------------------------------------------------
// run the diagnostic
bool FEResetTest::Run()
{
	FEBioModel* fem = dynamic_cast<FEBioModel*>(GetFEModel());

	// try to run the model
	if (fem->Solve() == false)
	{
		feLogEx(fem, "Failed to run model.");
		return false;
	}

	// reset the model
	if (fem->Reset() == false)
	{
		feLogEx(fem, "Failed to reset model.");
		return false;
	}

	// try to run it again
	if (fem->Solve() == false)
	{
		feLogEx(fem, "Failed to run model second time.");
		return false;
	}

	return true;
}
