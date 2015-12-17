#include "stdafx.h"
#include "FEBioModuleSection.h"

//-----------------------------------------------------------------------------
//! This function parses the Module section.
//! The Module defines the type of problem the user wants to solve (solid, heat, ...)
//! It is currently only used to allocate the FESolver.
void FEBioModuleSection::Parse(XMLTag &tag)
{
	int nversion = m_pim->Version();

	// For version 2.5 and up this tag can only be read in once and we 
	// determine that by inspecting the m_szmod member.
	if (nversion >= 0x0205)
	{
		const char* szmod = m_pim->m_szmod;
		if ((szmod != 0) && (strlen(szmod) > 0)) throw XMLReader::InvalidTag(tag);
	}

	// get the type attribute
	const char* szt = tag.AttributeValue("type");
	strcpy(m_pim->m_szmod, szt);
}
