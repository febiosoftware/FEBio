#include "stdafx.h"
#include "FEBioModuleSection.h"

//-----------------------------------------------------------------------------
//! This function parses the Module section.
//! The Module defines the type of problem the user wants to solve (solid, heat, ...)
//!
void FEBioModuleSection::Parse(XMLTag &tag)
{
	// get the type attribute
	const char* szt = tag.AttributeValue("type");

	if      (strcmp(szt, "solid"         ) == 0) m_pim->m_nstep_type = FE_SOLID;
	else if (strcmp(szt, "explicit-solid") == 0) m_pim->m_nstep_type = FE_EXPLICIT_SOLID;
	else if (strcmp(szt, "linear solid"  ) == 0) m_pim->m_nstep_type = FE_LINEAR_SOLID; 
	else if (strcmp(szt, "poro"          ) == 0) m_pim->m_nstep_type = FE_BIPHASIC;		// obsolete in 2.0
	else if (strcmp(szt, "biphasic"      ) == 0) m_pim->m_nstep_type = FE_BIPHASIC;
	else if (strcmp(szt, "solute"        ) == 0) m_pim->m_nstep_type = FE_POROSOLUTE;
	else if (strcmp(szt, "heat"          ) == 0) m_pim->m_nstep_type = FE_HEAT;
	else if (strcmp(szt, "heat-solid"    ) == 0) m_pim->m_nstep_type = FE_HEAT_SOLID;
	else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
}
