#include "stdafx.h"
#include "FEBioParametersSection.h"

//-----------------------------------------------------------------------------
//! Parse the Parameters section (new in version 2.0)
void FEBioParametersSection::Parse(XMLTag& tag)
{
	// make sure that the version is 2.x
	int nversion = m_pim->Version();
	if (nversion < 0x0200) throw XMLReader::InvalidTag(tag);

	// make sure there is something in here
	if (tag.isempty()) return;

	// get the xml reader
	XMLReader* pxml = tag.m_preader;

	// loop over parameters
	++tag;
	do
	{
		if (tag == "param")
		{
			const char* sz = tag.AttributeValue("name");
			m_pim->AddParameter(sz, tag.szvalue());
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
