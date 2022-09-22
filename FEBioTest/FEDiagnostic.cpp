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
#include "FEDiagnostic.h"
#include "FETangentDiagnostic.h"
#include "FEEASShellTangentDiagnostic.h"
#include "FEContactDiagnostic.h"
#include "FEPrintMatrixDiagnostic.h"
#include "FEPrintHBMatrixDiagnostic.h"
#include "FEMemoryDiagnostic.h"
#include "FEBiphasicTangentDiagnostic.h"
#include "FETiedBiphasicDiagnostic.h"
#include "FEMultiphasicTangentDiagnostic.h"
#include "FEFluidTangentDiagnostic.h"
#include "FEFluidFSITangentDiagnostic.h"
#include "FEPolarFluidTangentDiagnostic.h"
#include "FEContactDiagnosticBiphasic.h"
#include "FEMaterialTest.h"
#include "FECore/log.h"
#include "FEBioXML/FEBioControlSection.h"
#include "FEBioXML/FEBioMaterialSection.h"
#include "FEBioXML/FEBioGlobalsSection.h"
#include "FECore/FECoreKernel.h"
#include "FECore/FESolver.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEDiagnostic::FEDiagnostic(FEModel* fem) : FECoreClass(fem)
{

}

FEDiagnostic::~FEDiagnostic()
{

}

void FEDiagnostic::SetFileName(const std::string& fileName)
{
	m_file = fileName;
}

const std::string& FEDiagnostic::GetFileName()
{
	return m_file;
}

//-----------------------------------------------------------------------------
FEDiagnostic* FEDiagnosticImport::LoadFile(FEModel& fem, const char* szfile)
{
	m_pdia = 0;

	m_builder = new FEModelBuilder(fem);

	// Open the XML file
	XMLReader xml;
	if (xml.Open(szfile) == false) 
	{
		errf("FATAL ERROR: Failed opening input file %s\n\n", szfile);
		return 0;
	}

	// define file structure
	m_map.clear();
	m_map["Control" ] = new FEDiagnosticControlSection (this);
	m_map["Material"] = new FEBioMaterialSection       (this);
	m_map["Scenario"] = new FEDiagnosticScenarioSection(this);
    m_map["Globals" ] = new FEBioGlobalsSection        (this);

	FECoreKernel& fecore = FECoreKernel::GetInstance();

	// loop over all child tags
	try
	{
		// Find the root element
		XMLTag tag;
		if (xml.FindTag("febio_diagnostic", tag) == false) return 0;

		XMLAtt& att = tag.m_att[0];
        if      (att == "tangent test"            ) { fecore.SetActiveModule("solid"      ); m_pdia = new FETangentDiagnostic           (&fem); }
        else if (att == "shell tangent test"      ) { fecore.SetActiveModule("solid"      ); m_pdia = new FEEASShellTangentDiagnostic   (&fem); }
        else if (att == "contact test"            ) { fecore.SetActiveModule("solid"      ); m_pdia = new FEContactDiagnostic           (&fem); }
        else if (att == "print matrix"            ) { fecore.SetActiveModule("solid"      ); m_pdia = new FEPrintMatrixDiagnostic       (&fem); }
        else if (att == "print hbmatrix"          ) { fecore.SetActiveModule("solid"      ); m_pdia = new FEPrintHBMatrixDiagnostic     (&fem); }
        else if (att == "memory test"             ) { fecore.SetActiveModule("solid"      ); m_pdia = new FEMemoryDiagnostic            (&fem); }
        else if (att == "biphasic tangent test"   ) { fecore.SetActiveModule("biphasic"   ); m_pdia = new FEBiphasicTangentDiagnostic   (&fem); }
        else if (att == "biphasic contact test"   ) { fecore.SetActiveModule("biphasic"   ); m_pdia = new FEContactDiagnosticBiphasic   (&fem); }
        else if (att == "tied biphasic test"      ) { fecore.SetActiveModule("biphasic"   ); m_pdia = new FETiedBiphasicDiagnostic      (&fem); }
        else if (att == "multiphasic tangent test") { fecore.SetActiveModule("multiphasic"); m_pdia = new FEMultiphasicTangentDiagnostic(&fem); }
        else if (att == "fluid tangent test"      ) { fecore.SetActiveModule("fluid"      ); m_pdia = new FEFluidTangentDiagnostic      (&fem); }
        else if (att == "fluid-FSI tangent test"  ) { fecore.SetActiveModule("fluid-FSI"  ); m_pdia = new FEFluidFSITangentDiagnostic   (&fem); }
        else if (att == "polar fluid tangent test") { fecore.SetActiveModule("polar fluid"); m_pdia = new FEPolarFluidTangentDiagnostic (&fem); }
        else if (att == "material test"           ) { fecore.SetActiveModule("solid"      ); m_pdia = new FEMaterialTest                (&fem); }
		else
		{
			feLog("\nERROR: unknown diagnostic\n\n");
			return 0;
		}

        // keep a pointer to the fem object

		fem.SetCurrentStepIndex(0);
        
		// parse the file
		if (ParseFile(tag) == false) return nullptr;
	}
	catch (XMLReader::Error& e)
	{
		feLog("FATAL ERROR: %s\n", e.what());
		return 0;
	}
	catch (FEFileException& e)
	{
		feLog("FATAL ERROR: %s (line %d)\n", e.GetErrorString(), xml.GetCurrentLine());
		return 0;
	}
	catch (...)
	{
		feLog("FATAL ERROR: unrecoverable error (line %d)\n", xml.GetCurrentLine());
		return 0;
	}

	// close the XML file
	xml.Close();

	if (m_pdia) m_pdia->SetFileName(szfile);

	// we're done!
	return m_pdia;
}

//-----------------------------------------------------------------------------
void FEDiagnosticControlSection::Parse(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEAnalysis* pstep = fem.GetCurrentStep();

	++tag;
	do
	{
		if      (tag == "time_steps") tag.value(pstep->m_ntime);
		else if (tag == "step_size") { tag.value(pstep->m_dt0); fem.GetTime().timeIncrement = pstep->m_dt0; }
		else throw XMLReader::InvalidValue(tag);

		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEDiagnosticScenarioSection::Parse(XMLTag &tag)
{
	FEDiagnosticImport& dim = static_cast<FEDiagnosticImport&>(*GetFileReader());

	// get the diagnostic
	FEDiagnostic* pdia = dim.m_pdia;

	// find the type attribute
	XMLAtt& type = tag.Attribute("type");

	// create the scenario
	FEDiagnosticScenario* pscn = pdia->CreateScenario(type.cvalue());
	if (pscn == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", type.cvalue());

	// parse the parameter list
	FEParameterList& pl = pscn->GetParameterList();
	++tag;
	do
	{
		if (ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
