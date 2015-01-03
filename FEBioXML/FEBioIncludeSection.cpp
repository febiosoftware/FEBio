#include "stdafx.h"
#include "FEBioIncludeSection.h"

//-----------------------------------------------------------------------------
//! Parse the Include section (new in version 2.0)
//! This section includes the contents of another FEB file.
void FEBioIncludeSection::Parse(XMLTag& tag)
{
	// make sure that the version is 2.x
	int nversion = m_pim->Version();
	if (nversion < 0x0200) throw XMLReader::InvalidTag(tag);

	// see if we need to pre-pend a path
	char szin[512];
	strcpy(szin, tag.szvalue());
	char* ch = strrchr(szin, '\\');
	if (ch==0) ch = strrchr(szin, '/');
	if (ch==0)
	{
		// pre-pend the name with the input path
		sprintf(szin, "%s%s", m_pim->m_szpath, tag.szvalue());
	}

	// read the file
	if (m_pim->ReadFile(szin) == false)
		throw XMLReader::InvalidValue(tag);
}
