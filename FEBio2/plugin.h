#pragma once
#include <vector>

#ifdef WIN32
#include <Windows.h>
#undef RegisterClass
#undef GetFileTitle
#endif

#ifdef WIN32
typedef HMODULE FEBIO_PLUGIN_HANDLE;
#endif

#ifdef LINUX
typedef void* FEBIO_PLUGIN_HANDLE;
#endif

//-----------------------------------------------------------------------------
//! This class defines a FEBio plugin
class FEBioPlugin
{
public:
	FEBioPlugin();
	~FEBioPlugin();

	//! Try to load the library
	bool Load(const char* szfile);

	//! Unload the library
	void UnLoad();

private:
	FEBIO_PLUGIN_HANDLE		m_ph;
};

//-----------------------------------------------------------------------------
//! This class manages all the plugins
class FEBioPluginManager
{
public:
	//! Get the plugin manager
	static FEBioPluginManager* GetInstance();

	//! Load a plugin into memory
	bool LoadPlugin(const char* szfile);

	//! Clean up
	void DeleteThis();

private:
	std::vector<FEBioPlugin*>	m_Plugin;

private:
	FEBioPluginManager(){}
	~FEBioPluginManager();
	FEBioPluginManager(const FEBioPluginManager&){}
	static FEBioPluginManager*	m_pThis;
};
