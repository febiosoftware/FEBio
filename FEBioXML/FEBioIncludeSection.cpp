#include "stdafx.h"
#include "FEBioIncludeSection.h"

//-----------------------------------------------------------------------------
//! Parse the Include section (new in version 2.0)
//! This section includes the contents of another FEB file.
void FEBioIncludeSection::Parse(XMLTag& tag)
{
	// see if we need to pre-pend a path
	char szin[512];
	strcpy(szin, tag.szvalue());
	char* ch = strrchr(szin, '\\');
	if (ch==0) ch = strrchr(szin, '/');
	if (ch==0)
	{
		// pre-pend the name with the input path
		sprintf(szin, "%s%s", GetFileReader()->GetFilePath(), tag.szvalue());
	}

	// read the file
	if (GetFEBioImport()->ReadFile(szin, false) == false)
		throw XMLReader::InvalidValue(tag);
}
