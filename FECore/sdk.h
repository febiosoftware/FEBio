#pragma once
#include "FECoreKernel.h"

//------------------------------------------------------------------------------
// This macro should be used to start the definition of the plugin classes
#define BEGIN_PLUGIN_DEFINITION		static std::vector<FECoreFactory*>	_facs; \
FECORE_API void PluginInitialize(FECoreKernel& febio) { FECoreKernel::SetInstance(&febio);

//------------------------------------------------------------------------------
// Use this macro to define each new plugin class
#define REGISTER_PLUGIN(theClass, theID, theName)  \
	_facs.push_back(new FEPluginFactory_T<theClass, theID>(theName));

//------------------------------------------------------------------------------
// This macro ends the plugin definition section
#define END_PLUGIN_DEFINITION

#ifdef WIN32
	#define FECORE_PLUGIN extern "C" __declspec(dllexport)
#else
	#define FECORE_PLUGIN extern "C"
#endif
