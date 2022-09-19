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
#include "FEMaterialTest.h"
#include "FETangentDiagnostic.h"
#include <FECore/log.h>
#include <FECore/DataRecord.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FEDomain.h>

//-----------------------------------------------------------------------------
// Constructor
FEMaterialTest::FEMaterialTest(FEModel* fem) : FEDiagnostic(fem)
{
	m_szoutfile = "stress.txt";

	// make sure the correct module is active
	fem->SetActiveModule("solid");

	m_pscn = nullptr;
	m_strain = nullptr;
	m_stress = nullptr;

	// create an analysis step
	FEAnalysis* pstep = new FEAnalysis(fem);

	// create a new solver
	FESolver* pnew_solver = fecore_new<FESolver>("solid", fem);
	assert(pnew_solver);
	pstep->SetFESolver(pnew_solver);

	// keep a pointer to the fem object
	fem->AddStep(pstep);
	fem->SetCurrentStep(pstep);
}

//-----------------------------------------------------------------------------
FEMaterialTest::~FEMaterialTest()
{
	delete m_strain;
	delete m_stress;
}

//-----------------------------------------------------------------------------
void FEMaterialTest::SetOutputFileName(const char* szfilename)
{
	m_szoutfile = szfilename;
}

//-----------------------------------------------------------------------------
FEDiagnosticScenario* FEMaterialTest::CreateScenario(const std::string& sname)
{
	if (sname == "uni-axial"   ) return (m_pscn = new FETangentUniaxial(this));
	if (sname == "simple shear") return (m_pscn = new FETangentSimpleShear(this));
	return 0;
}

//-----------------------------------------------------------------------------
// Initialize the diagnostic. In this function we build the FE model depending
// on the scenario.
bool FEMaterialTest::Init()
{
	// make sure we have a scenario
	if (m_pscn == 0) return false;

	// initialize the scenario
	if (m_pscn->Init() == false) return false;

	FEModel* fem = GetFEModel();
	fem->AddCallback(cb, CB_INIT | CB_MAJOR_ITERS, this);

	m_strain = fecore_new<FELogElemData>("Ex", fem);
	m_stress = fecore_new<FELogElemData>("sx", fem);

	// add the stress output 
	if (m_szoutfile)
	{
		DataRecord* pdr = fecore_new<DataRecord>("element_data", fem);
		FEDomain& dom = fem->GetMesh().Domain(0);
		FEElementSet* es = new FEElementSet(fem);
		es->Create(&dom);
		std::vector<int> dummy;
		pdr->SetItemList(es, dummy);
		pdr->SetData("Ex;sx");
		pdr->SetComments(false);

		pdr->SetFileName(m_szoutfile);

		DataStore& ds = fem->GetDataStore();
		ds.AddRecord(pdr);
	}

	return FEDiagnostic::Init();
}

//-----------------------------------------------------------------------------
bool FEMaterialTest::cb()
{
	FEModel* fem = GetFEModel();
	FEDomain& dom= fem->GetMesh().Domain(0);
	FEElement& el = dom.ElementRef(0);
	
	double Ex = m_strain->value(el);
	double sx = m_stress->value(el);

	m_data.push_back(pair<double, double>(Ex, sx));

	return true;
}

//-----------------------------------------------------------------------------
// Run the tangent diagnostic. After we run the FE model, we calculate 
// the element stiffness matrix and compare that to a finite difference
// of the element residual.
bool FEMaterialTest::Run()
{
	// solve the problem
	FEModel& fem = *GetFEModel();
	fem.BlockLog();
	bool bret = fem.Solve();
	fem.UnBlockLog();
	if (bret == false)
	{
		feLogError("FEBio error terminated. Aborting diagnostic.\n");
		return false;
	}
	else
	{
		feLog("Material test completed.");
	}

	return true;
}
