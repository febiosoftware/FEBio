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
#include "FEDiagnostic.h"
#include <FECore/ElementDataRecord.h>

//-----------------------------------------------------------------------------
//! The material test diagnostics evaluates the stress-strain of a material. 
class FEMaterialTest : public FEDiagnostic
{
public:
	FEMaterialTest(FEModel* fem);
	~FEMaterialTest();

	FEDiagnosticScenario* CreateScenario(const std::string& sname);

	bool Init();

	bool Run();

	void SetOutputFileName(const char* szfilename);

	std::vector<std::pair<double, double> > GetOutputData() { return m_data; }

private:
	static bool cb(FEModel* fem, unsigned int when, void* pd) { return ((FEMaterialTest*)pd)->cb(); }
	bool cb();

public:
	FEDiagnosticScenario* m_pscn;
	
	const char* m_szoutfile;

	FELogElemData* m_strain;
	FELogElemData* m_stress;

	std::vector<std::pair<double, double> >	m_data;
};
