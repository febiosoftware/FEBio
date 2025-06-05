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
#include <FECore/DumpMemStream.h>

// This diagnostics tests the running and cold restart features
class FERestartDiagnostic : public FECoreTask
{
public:
	// constructor
	FERestartDiagnostic(FEModel* pfem);

	// initialize the diagnostic
	bool Init(const char* sz) override;

	// run the diagnostic
	bool Run() override;

public:
	bool	m_bok;
	bool	m_bfile;		// file or memory stream?
	char	m_szdmp[256];	// restart file name
	DumpMemStream	m_dmp;
};

// This diagnostics tests simple calls serialize after init
class FEQuickRestartDiagnostic : public FECoreTask
{
public:
	// constructor
	FEQuickRestartDiagnostic(FEModel* pfem);

	// initialize the diagnostic
	bool Init(const char* sz) override;

	// run the diagnostic
	bool Run() override;

public:
	char	m_szdmp[256];	// restart file name
};
