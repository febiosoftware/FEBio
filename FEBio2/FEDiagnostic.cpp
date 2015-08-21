// FEDiagnostic.cpp: implementation of the FEDiagnostic class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEDiagnostic.h"
#include "FETangentDiagnostic.h"
#include "FEContactDiagnostic.h"
#include "FEPrintMatrixDiagnostic.h"
#include "FEPrintHBMatrixDiagnostic.h"
#include "FEMemoryDiagnostic.h"
#include "FEBiphasicTangentDiagnostic.h"
#include "FECore/log.h"
#include "FEBioXML/FEBioControlSection.h"
#include "FEBioXML/FEBioMaterialSection.h"
#include "FEBioXML/FEBioGlobalsSection.h"
#include "FEBioMech/FESolidAnalysis.h"
#include "FECore/FECoreKernel.h"
#include "FECore/FESolver.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEDiagnostic::FEDiagnostic(FEModel& fem) : m_fem(fem)
{

}

FEDiagnostic::~FEDiagnostic()
{

}

//-----------------------------------------------------------------------------
FEDiagnostic* FEDiagnosticImport::LoadFile(FEModel& fem, const char* szfile)
{
	// Open the XML file
	XMLReader xml;
	if (xml.Open(szfile) == false) 
	{
		errf("FATAL ERROR: Failed opening input file %s\n\n", szfile);
		return 0;
	}

	// define file structure
	FEBioFileSectionMap map;
	map["Control" ] = new FEBioControlSection (this);
	map["Material"] = new FEBioMaterialSection(this);
	map["Scenario"] = new FEBioScenarioSection(this);
    map["Globals" ] = new FEBioGlobalsSection (this);

    FEAnalysis* pstep = 0;
    m_pfem = &fem;
    m_pdia = 0;
    
	// loop over all child tags
	try
	{
		// Find the root element
		XMLTag tag;
		if (xml.FindTag("febio_diagnostic", tag) == false) return 0;

		XMLAtt& att = tag.m_att[0];
        if      (att == "tangent test"  ) {
            m_pdia = new FETangentDiagnostic(fem);
            pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "solid", m_pfem);
        }
        else if (att == "contact test"  ) {
            m_pdia = new FEContactDiagnostic(fem);
            pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "solid", m_pfem);
        }
        else if (att == "print matrix"  ) {
            m_pdia = new FEPrintMatrixDiagnostic(fem);
            pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "solid", m_pfem);
        }
        else if (att == "print hbmatrix") {
            m_pdia = new FEPrintHBMatrixDiagnostic(fem);
            pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "solid", m_pfem);
        }
        else if (att == "memory test"   ) {
            m_pdia = new FEMemoryDiagnostic(fem);
            pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "solid", m_pfem);
        }
        else if (att == "biphasic tangent test") {
            m_pdia = new FEBiphasicTangentDiagnostic(fem);
            pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "biphasic", m_pfem);
            // create a new solver
            FESolver* pnew_solver = fecore_new<FESolver>(FESOLVER_ID, "biphasic", m_pfem);
            assert(pnew_solver);
            pnew_solver->m_bsymm = false;
            pstep->m_psolver = pnew_solver;
        }
		else
		{
			felog.printf("\nERROR: unknown diagnostic\n\n");
			return 0;
		}

        // keep a pointer to the fem object
        fem.AddStep(pstep);
        fem.m_nStep = 0;
        fem.SetCurrentStep(pstep);
        m_pStep = pstep;
        
		++tag;
		do
		{
			// parse the file
			FEBioFileSectionMap::iterator is = map.find(tag.Name());
			if (is != map.end()) is->second->Parse(tag);
			else throw XMLReader::InvalidTag(tag);

			// go to the next tag
			++tag;
		}
		while (!tag.isend());
	}
	catch (XMLReader::Error& e)
	{
		felog.printf("FATAL ERROR: %s (line %d)\n", e.GetErrorString(), xml.GetCurrentLine());
		return 0;
	}
	catch (FEBioImport::Exception& e)
	{
		felog.printf("FATAL ERROR: %s (line %d)\n", e.GetErrorString(), xml.GetCurrentLine());
		return 0;
	}
	catch (...)
	{
		felog.printf("FATAL ERROR: unrecoverable error (line %d)\n", xml.GetCurrentLine());
		return 0;
	}

	// close the XML file
	xml.Close();

	// we're done!
	return m_pdia;

}

//-----------------------------------------------------------------------------
void FEBioScenarioSection::Parse(XMLTag &tag)
{
	FEDiagnosticImport& dim = static_cast<FEDiagnosticImport&>(*m_pim);
	FETangentDiagnostic& td = static_cast<FETangentDiagnostic&>(*dim.m_pdia);
    FEBiphasicTangentDiagnostic& tbd = static_cast<FEBiphasicTangentDiagnostic&>(*dim.m_pdia);

	XMLAtt& type = tag.Attribute("type");
	if      (type == "uni-axial"   ) td.m_scn = FETangentDiagnostic::TDS_UNIAXIAL;
	else if (type == "simple shear") td.m_scn = FETangentDiagnostic::TDS_SIMPLE_SHEAR;
    else if (type == "biphasic uni-axial") tbd.m_scn = FEBiphasicTangentDiagnostic::TDS_BIPHASIC_UNIAXIAL;
	else throw XMLReader::InvalidAttributeValue(tag, "type", type.cvalue());
	++tag;
	do
	{
        if (tag == "strain") { tag.value(td.m_strain); }
        else if (tag == "solid_strain") { tag.value(tbd.m_strain); }
        else if (tag == "fluid_pressure") { tag.value(tbd.m_pressure); }
        else if (tag == "time_step") { tag.value(tbd.m_dt); }
		++tag;
	}
	while (!tag.isend());
}
