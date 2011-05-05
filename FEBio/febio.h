#pragma once

enum {
	FEBIO_MATERIAL
};

//-----------------------------------------------------------------------------
class IFEBioPlugin
{
public:
	virtual int GetPluginType() = 0;
	virtual void* CreateInstance() = 0;
};
