#include "stdafx.h"
#include "plugin.h"
#include "FECore/febio.h"
#include "log.h"

//-----------------------------------------------------------------------------
#ifndef WIN32
bool LoadPlugin(const char* szfile) { return false; }
#else
#include "windows.h"

typedef void (_cdecl *FEBIO_REGISTER_PLUGIN_FNC)(FEBioKernel&);

//-----------------------------------------------------------------------------
bool LoadPlugin(const char* szfile)
{
	FEBioKernel& febio = FEBioKernel::GetInstance();

	// load the library
	HMODULE hm = LoadLibraryA(szfile);
	if (hm == NULL) return false;

	// find the plugin's registration function
	FEBIO_REGISTER_PLUGIN_FNC pfnc = (FEBIO_REGISTER_PLUGIN_FNC) GetProcAddress(hm, "RegisterPlugin");
	if (pfnc == 0) return false;

	// allow the plugin to register itself with the framework
	pfnc(febio);

	// a-ok!
	return true;
}

#endif	// ifndef WIN32
