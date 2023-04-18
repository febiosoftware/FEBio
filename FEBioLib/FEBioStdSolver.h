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



#pragma once
#include <FECore/FECoreTask.h>
#include <FECore/FECoreKernel.h>

//-----------------------------------------------------------------------------
// This is the most commenly used task which will run a user-specified input 
// file. The results are stored in the logfile and the plotfile.
class FEBioStdSolver : public FECoreTask
{
public:
	FEBioStdSolver(FEModel* pfem);

	//! initialization
	bool Init(const char* szfile) override;

	//! Run the FE model
	bool Run() override;
};

//-----------------------------------------------------------------------------
// class for testing reverse communication interface of FEModel
class FEBioRCISolver : public FECoreTask
{
public:
	FEBioRCISolver(FEModel* fem);

	//! initialization
	bool Init(const char* szfile) override;

	//! Run the FE model
	bool Run() override;
};

//-----------------------------------------------------------------------------
// Configures the model for running in the nightly test suite. 
class FEBioTestSuiteTask : public FECoreTask
{
public:
	FEBioTestSuiteTask(FEModel* fem);

	//! initialization
	bool Init(const char* szfile) override;

	//! Run the FE model
	bool Run() override;
};
