#include "stdafx.h"
#include "plugin.h"
#include "FEBioLib/febio.h"
#include "fematerial.h"
#include "log.h"

//-----------------------------------------------------------------------------
#ifndef WIN32
bool LoadPlugin(const char* szfile) { return false; }
#else
#include "windows.h"

typedef void (_cdecl *FEBIO_REGISTER_PLUGIN_FNC)(FEBioKernel&);

extern FEBioKernel FEBio;

//-----------------------------------------------------------------------------
bool LoadPlugin(const char* szfile)
{
	// load the library
	HMODULE hm = LoadLibraryA(szfile);
	if (hm == NULL) return false;

	// find the plugin's registration function
	FEBIO_REGISTER_PLUGIN_FNC pfnc = (FEBIO_REGISTER_PLUGIN_FNC) GetProcAddress(hm, "RegisterPlugin");
	if (pfnc == 0) return false;

	// allow the plugin to register itself with the framework
	pfnc(FEBio);

	// a-ok!
	return true;
}

#endif	// ifndef WIN32
