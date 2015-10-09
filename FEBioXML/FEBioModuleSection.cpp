#include "stdafx.h"
#include "FEBioModuleSection.h"

//-----------------------------------------------------------------------------
//! This function parses the Module section.
//! The Module defines the type of problem the user wants to solve (solid, heat, ...)
//! It is currently only used to allocate the FESolver.
void FEBioModuleSection::Parse(XMLTag &tag)
{
	// get the type attribute
	const char* szt = tag.AttributeValue("type");
	strcpy(m_pim->m_szmod, szt);
}
