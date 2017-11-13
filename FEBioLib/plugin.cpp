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
	typedef unsigned int (*PLUGIN_GETSDKVERSION)();
	typedef void (*PLUGIN_INIT_FNC)(FECoreKernel&);
	typedef void (*PLUGIN_CLEANUP_FNC)();
	typedef int  (*PLUGIN_NUMCLASSES_FNC)();
	typedef FECoreFactory* (*PLUGIN_GETFACTORY_FNC)(int i);
	typedef void (*PLUGIN_VERSION_FNC)(int&,int&,int&);
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

	m_version.major = 0;
	m_version.minor = 0;
	m_version.patch = 0;
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
//! (2) Required plugin function PluginNumClasses not found (NOTE: as of 2.5 this error is not returned anymore since PluginNumClasses is obsolete)
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
	PLUGIN_GETSDKVERSION pf_sdk = (PLUGIN_GETSDKVERSION) FindPluginFunc(ph, "GetSDKVersion");
	if (pf_sdk == 0) return 5;

	// get the SDK version of the plugin
	unsigned int n = FE_SDK_VERSION;
	unsigned int version = pf_sdk();
	if (version != FE_SDK_VERSION) return 6;

	// find the numclasses function
	PLUGIN_NUMCLASSES_FNC pfnc_cnt = (PLUGIN_NUMCLASSES_FNC) FindPluginFunc(ph, "PluginNumClasses");
	
	// NOTE: As of 2.5, the PluginNumClasses is no longer required.
	//       If it is not defined, the PluginGetFactory will be called until null is returned.
	//       If it is defined, the behavior is as usual.
//	if (pfnc_cnt == 0) return 2;

	// find the GetFactory function
	PLUGIN_GETFACTORY_FNC pfnc_get = (PLUGIN_GETFACTORY_FNC) FindPluginFunc(ph, "PluginGetFactory");

	// NOTE: As of 2.6, the PluginGetFactory function is optional. This is because the 
	// REGISTER_FECORE_CLASS can be used to register a new plugin class in PluginInitialize.
//	if (pfnc_get == 0) return 3;
	
	// find the plugin's initialization function
	PLUGIN_INIT_FNC pfnc_init = (PLUGIN_INIT_FNC) FindPluginFunc(ph, "PluginInitialize");

	// find the optional plugin version function
	PLUGIN_VERSION_FNC pfnc_version = (PLUGIN_VERSION_FNC) FindPluginFunc(ph, "GetPluginVersion");
	if (pfnc_version) pfnc_version(m_version.major, m_version.minor, m_version.patch);

	// call the (optional) initialization function
	FECoreKernel& febio = FECoreKernel::GetInstance();
	if (pfnc_init) 
	{
		pfnc_init(febio);
		febio.SetActiveModule(0);
	}

	// find out how many classes there are in this plugin
	// This is only called when the PluginGetFactory function was found. 
	// (As of 2.6, this is no longer a required function)
	if (pfnc_get)
	{
		if (pfnc_cnt)
		{
			int NC = pfnc_cnt();
			// call the get factory functions
			for (int i=0; i<NC; ++i)
			{
				FECoreFactory* pfac = pfnc_get(i);
				if (pfac) febio.RegisterClass(pfac);
			}
		}
		else
		{
			// As of 2.5, the PluginNumClasses is no longer required. 
			// In this case, the PluginGetFactory is called until null is returned.
			FECoreFactory* pfac = 0;
			int i = 0;
			do
			{
				pfac = pfnc_get(i);
				if (pfac)
				{
					febio.RegisterClass(pfac);
					i++;
				}
			}
			while (pfac);
		}
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
		PLUGIN_CLEANUP_FNC pfnc = (PLUGIN_CLEANUP_FNC) FindPluginFunc(m_ph, "PluginCleanup");
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
int FEBioPluginManager::LoadPlugin(const char* szfile, PLUGIN_INFO& info)
{
	// create a new plugin object
	FEBioPlugin* pdll = new FEBioPlugin;

	info.bloaded = false;
	info.major = 0;
	info.minor = 0;
	info.patch = 0;

	// try to load the plugin
	int nerr = pdll->Load(szfile);

	// add it to the list or delete it if error
	if (nerr == 0) 
	{
		FEBioPlugin::Version version = pdll->GetVersion();

		info.bloaded = true;
		info.major = version.major;
		info.minor = version.minor;
		info.patch = version.patch;

		m_Plugin.push_back(pdll);
	}
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
