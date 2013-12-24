// dllmain.cpp : Defines the entry point for the DLL application.
#include "FECore/FECoreKernel.h"
#include "FECore/FECoreFactory.h"
#include "FECore/Logfile.h"
#include "FEBioOpt.h"

#ifdef WIN32
	#define DLL_EXPORT __declspec(dllexport)
#else
	#define DLL_EXPORT
#endif

#ifdef WIN32
#include "targetver.h"
#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
#include <windows.h>

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
	switch (ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
	return TRUE;
}

#ifdef RegisterClass
#undef RegisterClass
#endif

#endif // WIN32

//-----------------------------------------------------------------------------
// the optimization task factory
class FEBioOptFactory : public FECoreFactory
{
public:
	FEBioOptFactory() : FECoreFactory(FETASK_ID, "optimize"){}
	void* Create(FEModel* pfem) { return new FEBioOpt(pfem); }
};

FEBioOptFactory	febioopt_factory;

//-----------------------------------------------------------------------------
// keep a copy of the FECoreKernel
FECoreKernel*	pFEBio;

//-----------------------------------------------------------------------------
extern "C" DLL_EXPORT void PluginInitialize(FECoreKernel& febio)
{
	pFEBio = &febio;
}

//-----------------------------------------------------------------------------
extern "C" DLL_EXPORT int PluginNumClasses()
{
	return 1;
}

//-----------------------------------------------------------------------------
extern "C" DLL_EXPORT FECoreFactory* PluginGetFactory(int i)
{
	if (i==0) return &febioopt_factory;
	return 0;
}

//-----------------------------------------------------------------------------
extern "C" DLL_EXPORT void PluginCleanup()
{

}
