#include "stdafx.h"
#include "plugin.h"
#include "FECore/FECoreKernel.h"
#include "FECore/log.h"

#ifdef WIN32
#include <direct.h>
#endif
#ifdef LINUX
#include <dlfcn.h>
#endif
#ifdef __APPLE__
#include <dlfcn.h>
#endif

//=============================================================================
// Typedefs of the plugin functions.
// These are the functions that the plugin must implement
extern "C" {
	typedef unsigned int (*FEPLUGIN_GETSDKVERSION)();
	typedef void (*FEPLUGIN_INIT_FNC)(FECoreKernel&);
	typedef void (*FEPLUGIN_CLEANUP_FNC)();
	typedef int  (*FEPLUGIN_NUMCLASSES_FNC)();
	typedef FECoreFactory* (*FEPLUGIN_GETFACTORY_FNC)(int i);
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
#ifdef __APPLE__
FEBIO_PLUGIN_HANDLE LoadPlugin(const char* szfile) { return dlopen(szfile, RTLD_NOW); }
void* FindPluginFunc(FEBIO_PLUGIN_HANDLE ph, const char* szfunc) { return dlsym(ph, szfunc); }
#endif

//=============================================================================
// FEBioPlugin
//=============================================================================
FEBioPlugin::FEBioPlugin()
{
	m_ph = 0;
	m_szname[0] = 0;
}

//-----------------------------------------------------------------------------
FEBioPlugin::~FEBioPlugin()
{
	if (m_ph) UnLoad();
}

//-----------------------------------------------------------------------------
void FEBioPlugin::SetNameFromFilePath(const char* szfile)
{
	const char* ch = strrchr(szfile, '\\');
	if (ch==0) 
	{
		ch = strrchr(szfile, '/'); 
		if (ch==0) ch = szfile; else ch++;
	} else ch++;
	if (ch) strcpy(m_szname, ch);
}

//-----------------------------------------------------------------------------
//! This function tries to load a plugin from file.
//! \return Return values are: 
//! (0) success, 
//! (1) Failed to load the file, 
//! (2) Required plugin function PluginNumClasses not found,
//! (3) Required plugin function PluginGetFactory not found,
//! (4) Invalid number of classes returned by PluginNumClasses.
//! (5) Required plugin function GetSDKVersion not found.
//! (6) Invalid SDK version.
int FEBioPlugin::Load(const char* szfile)
{
	// Make sure the plugin is not loaded already
	assert(m_ph == 0);
	if (m_ph) return 0;

	// set the file name as the plugin name
	SetNameFromFilePath(szfile);

	// load the library
	FEBIO_PLUGIN_HANDLE ph = LoadPlugin(szfile);
	if (ph == NULL) return 1;

	// get the GetSDKVersion function
	FEPLUGIN_GETSDKVERSION pf_sdk = (FEPLUGIN_GETSDKVERSION) FindPluginFunc(ph, "GetSDKVersion");
	if (pf_sdk == 0) return 5;

	// get the SDK version of the plugin
	unsigned int n = FE_SDK_VERSION;
	unsigned int version = pf_sdk();
	if (version != FE_SDK_VERSION) return 6;

	// find the numclasses function
	FEPLUGIN_NUMCLASSES_FNC pfnc_cnt = (FEPLUGIN_NUMCLASSES_FNC) FindPluginFunc(ph, "PluginNumClasses");
	if (pfnc_cnt == 0) return 2;

	// find the GetFactory function
	FEPLUGIN_GETFACTORY_FNC pfnc_get = (FEPLUGIN_GETFACTORY_FNC) FindPluginFunc(ph, "PluginGetFactory");
	if (pfnc_get == 0) return 3;
	
	// find the plugin's initialization function
	FEPLUGIN_INIT_FNC pfnc_init = (FEPLUGIN_INIT_FNC) FindPluginFunc(ph, "PluginInitialize");

	// call the (optional) initialization function
	FECoreKernel& febio = FECoreKernel::GetInstance();
	if (pfnc_init) pfnc_init(febio);

	// find out how many classes there are in this plugin
	int NC = pfnc_cnt();
	if (NC < 0) return 4;

	// call the get factory functions
	for (int i=0; i<NC; ++i)
	{
		FECoreFactory* pfac = pfnc_get(i);
		if (pfac) febio.RegisterClass(pfac);
	}

	// If we get here everything seems okay so let's store the handle
	m_ph = ph;

	// a-ok!
	return 0;
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
int FEBioPluginManager::Plugins()
{ 
	return (int) m_Plugin.size(); 
}

//-----------------------------------------------------------------------------
const FEBioPlugin& FEBioPluginManager::GetPlugin(int i)
{
	return *(m_Plugin[i]);
}

//-----------------------------------------------------------------------------
//! This function tries to load a plugin and returns the error code from the
//! FEBioPlugin::Load function. 
//! \return Returns zero on success, nonzero on failure.
//! \sa FEBioPlugin::Load
int FEBioPluginManager::LoadPlugin(const char* szfile)
{
	// create a new plugin object
	FEBioPlugin* pdll = new FEBioPlugin;

	// try to load the plugin
	int nerr = pdll->Load(szfile);

	// add it to the list or delete it if error
	if (nerr == 0) m_Plugin.push_back(pdll);
	else delete pdll;

	// pass error code to caller
	return nerr;
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
