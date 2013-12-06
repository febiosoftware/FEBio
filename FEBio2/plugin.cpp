#include "stdafx.h"
#include "plugin.h"
#include "FECore/febio.h"
#include "FECore/log.h"

#ifdef WIN32
#include <direct.h>
#endif
#ifdef LINUX
#include <dlfcn.h>
#endif

//=============================================================================
// Typedefs of the plugin functions.
// These are the functions that the plugin must implement
extern "C" {
	typedef void (*FEPLUGIN_INIT_FNC)(FEBioKernel&);
	typedef void (*FEPLUGIN_CLEANUP_FNC)();
	typedef int  (*FEPLUGIN_NUMCLASSES_FNC)();
	typedef FEBioFactory* (*FEPLUGIN_GETFACTORY_FNC)(int i);
}

//=============================================================================
// Wrappers for system calls
#ifdef WIN32
FEBIO_PLUGIN_HANDLE LoadPlugin(const char* szfile) { return LoadLibraryA(szfile); }
void* FindPluginFunc(FEBIO_PLUGIN_HANDLE ph, const char* szfunc) { return GetProcAddress(ph, szfunc); }
#endif
#ifdef LINUX
FEBIO_PLUGIN_HANDLE LoadPlugin(const char* szfile) { return dlopen(szfile, RTLD_NOW); }
void* FindPluginFunc(FEBIO_PLUGIN_HANDLE ph, const char* szfunc) { return dlsym(ph, szfunc); }
#endif

//=============================================================================
// FEBioPlugin
//=============================================================================
FEBioPlugin::FEBioPlugin()
{
	m_ph = 0;
}

//-----------------------------------------------------------------------------
FEBioPlugin::~FEBioPlugin()
{
	if (m_ph) UnLoad();
}

//-----------------------------------------------------------------------------
bool FEBioPlugin::Load(const char* szfile)
{
	// Make sure the plugin is not loaded already
	assert(m_ph == 0);
	if (m_ph) return true;

	// load the library
	FEBIO_PLUGIN_HANDLE ph = LoadPlugin(szfile);
	if (ph == NULL) return false;

	// find the numclasses function
	FEPLUGIN_NUMCLASSES_FNC pfnc_cnt = (FEPLUGIN_NUMCLASSES_FNC) FindPluginFunc(ph, "PluginNumClasses");
	if (pfnc_cnt == 0) return false;

	// find the GetFactory function
	FEPLUGIN_GETFACTORY_FNC pfnc_get = (FEPLUGIN_GETFACTORY_FNC) FindPluginFunc(ph, "PluginGetFactory");
	if (pfnc_get == 0) return false;
	
	// find the plugin's initialization function
	FEPLUGIN_INIT_FNC pfnc_init = (FEPLUGIN_INIT_FNC) FindPluginFunc(ph, "PluginInitialize");

	// call the (optional) initialization function
	FEBioKernel& febio = FEBioKernel::GetInstance();
	if (pfnc_init) pfnc_init(febio);

	// find out how many classes there are in this plugin
	int NC = pfnc_cnt();
	if (NC < 0) return false;

	// call the get factory functions
	for (int i=0; i<NC; ++i)
	{
		FEBioFactory* pfac = pfnc_get(i);
		if (pfac) febio.RegisterClass(pfac);
	}

	// If we get here everything seems okay so let's store the handle
	m_ph = ph;

	// a-ok!
	return true;

}

//-----------------------------------------------------------------------------
void FEBioPlugin::UnLoad()
{
	if (m_ph)
	{
		// find the plugin's cleanup function
		FEPLUGIN_CLEANUP_FNC pfnc = (FEPLUGIN_CLEANUP_FNC) FindPluginFunc(m_ph, "PluginCleanup");
		if (pfnc) pfnc();

		// TODO: should I figure out how to actually unload the plugin from memory?
		m_ph = 0;
	}
}

//=============================================================================
// FEBioPluginManager
//=============================================================================

FEBioPluginManager* FEBioPluginManager::m_pThis = 0;

FEBioPluginManager* FEBioPluginManager::GetInstance()
{
	if (m_pThis == 0) m_pThis = new FEBioPluginManager;
	return m_pThis;
}

//-----------------------------------------------------------------------------
FEBioPluginManager::~FEBioPluginManager()
{
	for (size_t i = 0; i < m_Plugin.size(); ++i) delete m_Plugin[i];
	m_Plugin.clear();
}

//-----------------------------------------------------------------------------
void FEBioPluginManager::DeleteThis()
{
	delete m_pThis;
	m_pThis = 0;
}

//-----------------------------------------------------------------------------
bool FEBioPluginManager::LoadPlugin(const char* szfile)
{
	// create a new plugin object
	FEBioPlugin* pdll = new FEBioPlugin;

	// try to load the plugin
	if (pdll->Load(szfile) == false)
	{
		delete pdll;
		return false;
	}
	else
	{
		m_Plugin.push_back(pdll);
		return true;
	}
}

/*
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
*/
