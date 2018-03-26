// FEDiagnostic.cpp: implementation of the FEDiagnostic class.
//
//////////////////////////////////////////////////////////////////////

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
#include "FEContactDiagnosticBiphasic.h"
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
	m_pdia = 0;

	if (Load(fem, szfile) == false) return 0;

	return m_pdia;
}

//-----------------------------------------------------------------------------
bool FEDiagnosticImport::Parse(const char* szfile)
{
	// Open the XML file
	XMLReader xml;
	if (xml.Open(szfile) == false) 
	{
		errf("FATAL ERROR: Failed opening input file %s\n\n", szfile);
		return 0;
	}

	// define file structure
	FEFileSectionMap map;
	map["Control" ] = new FEDiagnosticControlSection (this);
	map["Material"] = new FEBioMaterialSection       (this);
	map["Scenario"] = new FEDiagnosticScenarioSection(this);
    map["Globals" ] = new FEBioGlobalsSection        (this);

	FEModel& fem = *GetFEModel();

	// loop over all child tags
	try
	{
		// Find the root element
		XMLTag tag;
		if (xml.FindTag("febio_diagnostic", tag) == false) return 0;

		XMLAtt& att = tag.m_att[0];
        if      (att == "tangent test"            ) m_pdia = new FETangentDiagnostic           (fem);
        else if (att == "shell tangent test"      ) m_pdia = new FEEASShellTangentDiagnostic   (fem);
        else if (att == "contact test"            ) m_pdia = new FEContactDiagnostic           (fem);
        else if (att == "print matrix"            ) m_pdia = new FEPrintMatrixDiagnostic       (fem);
        else if (att == "print hbmatrix"          ) m_pdia = new FEPrintHBMatrixDiagnostic     (fem);
        else if (att == "memory test"             ) m_pdia = new FEMemoryDiagnostic            (fem);
        else if (att == "biphasic tangent test"   ) m_pdia = new FEBiphasicTangentDiagnostic   (fem);
        else if (att == "biphasic contact test"   ) m_pdia = new FEContactDiagnosticBiphasic   (fem);
        else if (att == "tied biphasic test"      ) m_pdia = new FETiedBiphasicDiagnostic      (fem);
        else if (att == "multiphasic tangent test") m_pdia = new FEMultiphasicTangentDiagnostic(fem);
        else if (att == "fluid tangent test"      ) m_pdia = new FEFluidTangentDiagnostic      (fem);
        else if (att == "fluid-FSI tangent test"  ) m_pdia = new FEFluidFSITangentDiagnostic   (fem);
		else
		{
			felog.printf("\nERROR: unknown diagnostic\n\n");
			return 0;
		}

        // keep a pointer to the fem object

		fem.SetCurrentStepIndex(0);
        
		// parse the file
		map.Parse(tag);
	}
	catch (XMLReader::Error& e)
	{
		felog.printf("FATAL ERROR: %s (line %d)\n", e.GetErrorString(), xml.GetCurrentLine());
		return 0;
	}
	catch (FEFileException& e)
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
	return true;
}

//-----------------------------------------------------------------------------
void FEDiagnosticControlSection::Parse(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEAnalysis* pstep = new FEAnalysis(&fem);

	++tag;
	do
	{
		if      (tag == "time_steps") tag.value(pstep->m_ntime);
		else if (tag == "step_size" ) { tag.value(pstep->m_dt0); pstep->m_dt = pstep->m_dt0; }
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
