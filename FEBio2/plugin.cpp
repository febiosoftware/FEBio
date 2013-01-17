#include "stdafx.h"
#include "plugin.h"
#include "FECore/febio.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
#ifdef WIN32

#include "windows.h"
#include <direct.h>

typedef void (_cdecl *FEBIO_REGISTER_PLUGIN_FNC)(FEBioKernel&);

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

//-----------------------------------------------------------------------------
// Import all the plugins from a folder. The szdir is actually a folder name
// plus a wildcard file reference (e.g. C:\folder\*.dll)
bool LoadPluginFolder(const char* szdir)
{
	WIN32_FIND_DATAA FileData;
	HANDLE hFind = FindFirstFileA(szdir, &FileData);

	LoadPlugin(FileData.cFileName);
	printf("Plugin \"%s\" loaded successfully\n", FileData.cFileName);

	while (FindNextFileA(hFind, &FileData)) 
	{
		LoadPlugin(FileData.cFileName);
		printf("Plugin \"%s\" loaded successfully\n", FileData.cFileName);
	}

	return true;
}

#endif	// ifdef WIN32


//-----------------------------------------------------------------------------
#ifdef LINUX
#include <dlfcn.h>

extern "C" {
typedef void (*FEBIO_REGISTER_PLUGIN_FNC)(FEBioKernel&);
}

bool LoadPlugin(const char* szfile)
{
  // load the library
  void* hlib = dlopen(szfile, RTLD_LAZY);
  if (hlib == NULL) return false;

  // find the plugin's registration function
  FEBIO_REGISTER_PLUGIN_FNC pfnc = (FEBIO_REGISTER_PLUGIN_FNC) dlsym(hlib, "RegisterPlugin");
  if (pfnc == NULL) return false;

  // allow the plugin to register itself with the framework
  FEBioKernel& febio = FEBioKernel::GetInstance();
  pfnc(febio);

  return true;
}

bool LoadPluginFolder(const char* szdir)
{
	return false;
}

#endif // ifdef LINUX
