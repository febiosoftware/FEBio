#include "stdafx.h"
#include "plugin.h"
#include "febio.h"
#include "fematerial.h"
#include "log.h"

//-----------------------------------------------------------------------------
#ifndef WIN32
bool LoadPlugin(const char* szfile) { return false; }
#else
#include "windows.h"

typedef IFEBioPlugin* (_cdecl *FEBIO_PLUGIN_FNC)();

//-----------------------------------------------------------------------------
bool LoadPlugin(const char* szfile)
{
	HMODULE hm = LoadLibraryA(szfile);
	if (hm == NULL) return false;

	FEBIO_PLUGIN_FNC pfnc = (FEBIO_PLUGIN_FNC) GetProcAddress(hm, "CreatePlugin");
	if (pfnc == 0) return false;

	// get the plugin's class descriptor
	IFEBioPlugin* pcd = pfnc();
	if (pcd == 0) return false;

	switch (pcd->GetPluginType())
	{
	case FEBIO_MATERIAL:
		{
			clog.printf("Inside switch\n");
		};
		break;
	default:
		return false;
	}

	// a-ok!
	return true;
}

#endif	// ifndef WIN32
