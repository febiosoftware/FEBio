#include "stdafx.h"
#include "FEBioStepSection.h"
#include "FEBioModuleSection.h"
#include "FEBioControlSection.h"
#include "FEBioConstraintsSection.h"
#include "FEBioBoundarySection.h"
#include "FEBioLoadsSection.h"
#include "FEBioContactSection.h"
#include "FEBioInitialSection.h"

//=============================================================================
//
//                         S T E P   S E C T I O N
//
//=============================================================================

void FEBioStepSection::Parse(XMLTag& tag)
{
	// reset the step pointer
	if (m_pim->m_nsteps != 0) m_pim->m_pStep = 0;

	// increase the step section counter
	++m_pim->m_nsteps;

	FEBioFileSectionMap Map;
	if (m_pim->Version() < 0x0205)
	{
		// Module section must be defined at the top of the file
		// for version 2.5 and cannot be redefined
		Map["Module"     ] = new FEBioModuleSection     (m_pim);
	}
	Map["Control"    ] = new FEBioControlSection    (m_pim);
	Map["Constraints"] = new FEBioConstraintsSection(m_pim);
	Map["Boundary"   ] = new FEBioBoundarySection   (m_pim);
	Map["Loads"      ] = new FEBioLoadsSection      (m_pim);
	Map["Initial"    ] = new FEBioInitialSection    (m_pim);

	if (m_pim->Version() >= 0x0200)
	{
		Map["Contact"] = new FEBioContactSection(m_pim);
	}

	++tag;
	do
	{
		std::map<string, FEBioFileSection*>::iterator is = Map.find(tag.Name());
		if (is != Map.end()) is->second->Parse(tag);
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}
