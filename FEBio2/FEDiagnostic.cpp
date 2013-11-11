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
#include "FECore/log.h"
#include "FEBioXML/FEBioControlSection.h"
#include "FEBioXML/FEBioMaterialSection.h"
#include "FEBioMech/FESolidAnalysis.h"

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

	// keep a pointer to the fem object
	m_pfem = &fem;
	m_pdia = 0;

	FEAnalysis* pstep = new FESolidAnalysis(m_pfem);
	fem.AddStep(pstep);
	fem.m_nStep = 0;
	fem.SetCurrentStep(pstep);
	m_pStep = pstep;

	// define file structure
	FEBioFileSectionMap map;
	map["Control" ] = new FEBioControlSection (this);
	map["Material"] = new FEBioMaterialSection(this);
	map["Scenario"] = new FEBioScenarioSection(this);

	// loop over all child tags
	try
	{
		// Find the root element
		XMLTag tag;
		if (xml.FindTag("febio_diagnostic", tag) == false) return false;

		XMLAtt& att = tag.m_att[0];
		if      (att == "tangent test"  ) m_pdia = new FETangentDiagnostic(fem);
		else if (att == "contact test"  ) m_pdia = new FEContactDiagnostic(fem);
		else if (att == "print matrix"  ) m_pdia = new FEPrintMatrixDiagnostic(fem);
		else if (att == "print hbmatrix") m_pdia = new FEPrintHBMatrixDiagnostic(fem);
		else if (att == "memory test"   ) m_pdia = new FEMemoryDiagnostic(fem);
		else 
		{
			clog.printf("\nERROR: unknown diagnostic\n\n");
			return 0;
		}

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
	catch (InvalidVersion)
	{
		clog.printbox("FATAL ERROR", "Invalid version for FEBio specification.");
		return false;
	}
	catch (InvalidMaterial e)
	{
		clog.printbox("FATAL ERROR:", "Element %d has an invalid material type.", e.m_nel);
		return false;
	}
	catch (XMLReader::XMLSyntaxError)
	{
		clog.printf("FATAL ERROR: Syntax error (line %d)\n", xml.GetCurrentLine());
		return false;
	}
	catch (XMLReader::InvalidTag e)
	{
		clog.printf("FATAL ERROR: unrecognized tag \"%s\" (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::InvalidAttributeValue e)
	{
		const char* szt = e.tag.m_sztag;
		const char* sza = e.szatt;
		const char* szv = e.szval;
		int l = e.tag.m_nstart_line;
		clog.printf("FATAL ERROR: unrecognized value \"%s\" for attribute \"%s.%s\" (line %d)\n", szv, szt, sza, l);
		return false;
	}
	catch (XMLReader::InvalidValue e)
	{
		clog.printf("FATAL ERROR: the value for tag \"%s\" is invalid (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::MissingAttribute e)
	{
		clog.printf("FATAL ERROR: Missing attribute \"%s\" of tag \"%s\" (line %d)\n", e.szatt, e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::UnmatchedEndTag e)
	{
		const char* sz = e.tag.m_szroot[e.tag.m_nlevel];
		clog.printf("FATAL ERROR: Unmatched end tag for \"%s\" (line %d)\n", sz, e.tag.m_nstart_line);
		return false;
	}
	catch (...)
	{
		clog.printf("FATAL ERROR: unrecoverable error (line %d)\n", xml.GetCurrentLine());
		return false;
	}

	// close the XML file
	xml.Close();

	// we're done!
	return m_pdia;

}

//-----------------------------------------------------------------------------
void FEBioScenarioSection::Parse(XMLTag &tag)
{
	FEDiagnosticImport& dim = dynamic_cast<FEDiagnosticImport&>(*m_pim);
	FETangentDiagnostic& td = dynamic_cast<FETangentDiagnostic&>(*dim.m_pdia);

	XMLAtt& type = tag.Attribute("type");
	if      (type == "uni-axial"   ) td.m_scn = FETangentDiagnostic::TDS_UNIAXIAL;
	else if (type == "simple shear") td.m_scn = FETangentDiagnostic::TDS_SIMPLE_SHEAR;
	else throw XMLReader::InvalidAttributeValue(tag, "type", type.cvalue());
	++tag;
	do
	{
		if (tag == "strain") tag.value(td.m_strain);
		++tag;
	}
	while (!tag.isend());
}
