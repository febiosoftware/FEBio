#include "stdafx.h"
#include "FEBioParametersSection.h"

//-----------------------------------------------------------------------------
//! Parse the Parameters section
void FEBioParametersSection::Parse(XMLTag& tag)
{
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
			GetFileReader()->AddFileParameter(sz, tag.szvalue());
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
