#pragma once
#include <vector>
#include <string>

#ifdef WIN32
#include <Windows.h>
#undef RegisterClass
#undef GetFileTitle
#endif

#include "febiolib_api.h"

#ifdef WIN32
typedef HMODULE FEBIO_PLUGIN_HANDLE;
#endif

#ifdef LINUX
typedef void* FEBIO_PLUGIN_HANDLE;
#endif

#ifdef __APPLE__
typedef void* FEBIO_PLUGIN_HANDLE;
#endif

class FECoreFactory;

//-----------------------------------------------------------------------------
struct PLUGIN_INFO
{
	bool	bloaded;	// loaded successfully
	int		major;	// major version number
	int		minor;	// minor version number
	int		patch;	// patch versin number
};

//-----------------------------------------------------------------------------
//! This class defines a FEBio plugin
class FEBIOLIB_API FEBioPlugin
{
public:
	struct Version
	{
		int major;	// major version number
		int	minor;	// minor version number
		int	patch;	// patch versin number
	};

public:
	FEBioPlugin();
	~FEBioPlugin();

	//! Try to load the library
	int Load(const char* szfile);

	//! Unload the library
	void UnLoad();

	//! get the plugin name
	const char* GetName() const { return m_szname; }

	//! return the version info
	Version GetVersion() const { return m_version; }

	//! get the plugin's path
	std::string GetFilePath() const { return m_filepath; }

protected:
	void SetNameFromFilePath(const char* szfile);

private:
	std::string					m_filepath;
	char						m_szname[1024];
	int							m_allocater_id;
	Version						m_version;
	FEBIO_PLUGIN_HANDLE			m_ph;
};

//-----------------------------------------------------------------------------
//! This class manages all the plugins
class FEBIOLIB_API FEBioPluginManager
{
public:
	//! Get the plugin manager
	static FEBioPluginManager* GetInstance();

	//! Load a plugin into memory
	int LoadPlugin(const char* szfile, PLUGIN_INFO& info);

	//! unload a plugin from memory
	bool UnloadPlugin(int n);
	bool UnloadPlugin(const std::string& name);
	void UnloadAllPlugins();

	//! Clean up
	void DeleteThis();

	//! return the number of plugins loaded
	int Plugins();

	//! return an instance of a plugin
	const FEBioPlugin& GetPlugin(int i);

private:
	std::vector<FEBioPlugin*>	m_Plugin;

private:
	FEBioPluginManager(){}
	~FEBioPluginManager();
	FEBioPluginManager(const FEBioPluginManager&){}
	static FEBioPluginManager*	m_pThis;
};
