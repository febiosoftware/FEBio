#include "stdafx.h"
#include "FECore.h"

//-----------------------------------------------------------------------------
//! The FECore namespace encapsulates all classes that belong to the FECore library

//! The FECore library is a collection of tools that simplify the developement
//! of FE software. 

namespace FECore
{

#define FECORE_VERSION		0
#define FECORE_SUBVERSION	1

static char fecore_str[4] = {'0'+FECORE_VERSION, '.', '0'+FECORE_SUBVERSION };

void get_version(int& version, int& subversion)
{
	version = FECORE_VERSION;
	subversion = FECORE_SUBVERSION;
}

const char* get_version_string()
{
	return fecore_str;
}

} // namespace FECore
