/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



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

	//! Get the allocator ID
	int GetAllocatorID() const { return m_allocater_id; }

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
