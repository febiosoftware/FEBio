// FEDiagnostic.cpp: implementation of the FEDiagnostic class.
//
//////////////////////////////////////////////////////////////////////

#include "FEDiagnostic.h"
#include "FETangentDiagnostic.h"
#include "FEContactDiagnostic.h"
#include "FEPrintMatrixDiagnostic.h"
#include "FEPrintHBMatrixDiagnostic.h"
#include "FEMemoryDiagnostic.h"
#include "FEBiphasicTangentDiagnostic.h"
#include "FEMultiphasicTangentDiagnostic.h"
#include "FECore/log.h"
#include "FEBioXML/FEBioControlSection.h"
#include "FEBioXML/FEBioMaterialSection.h"
#include "FEBioXML/FEBioGlobalsSection.h"
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
        if      (att == "tangent test"            ) m_pdia = new FETangentDiagnostic           (fem);
        else if (att == "contact test"            ) m_pdia = new FEContactDiagnostic           (fem);
        else if (att == "print matrix"            ) m_pdia = new FEPrintMatrixDiagnostic       (fem);
        else if (att == "print hbmatrix"          ) m_pdia = new FEPrintHBMatrixDiagnostic     (fem);
        else if (att == "memory test"             ) m_pdia = new FEMemoryDiagnostic            (fem);
        else if (att == "biphasic tangent test"   ) m_pdia = new FEBiphasicTangentDiagnostic   (fem);
        else if (att == "multiphasic tangent test") m_pdia = new FEMultiphasicTangentDiagnostic(fem);
		else
		{
			felog.printf("\nERROR: unknown diagnostic\n\n");
			return 0;
		}

        // keep a pointer to the fem object
        m_pStep = fem.GetCurrentStep();
        
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

	// get the diagnostic
	FEDiagnostic* pdia = dim.m_pdia;

	// find the type attribute
	XMLAtt& type = tag.Attribute("type");

	// create the scenario
	FEDiagnosticScenario* pscn = pdia->CreateScenario(type.cvalue());

	// parse the parameter list
	FEParameterList& pl = pscn->GetParameterList();
	++tag;
	do
	{
		if (m_pim->ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
