// dllmain.cpp : Defines the entry point for the DLL application.
#include "FECore/febio.h"
#include "FECore/FEBioFactory.h"
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
class FEBioOptFactory : public FEBioFactory_T<FEBioTask>
{
public:
	FEBioOptFactory() : FEBioFactory_T<FEBioTask>("optimize"){}
	bool IsType(FEBioTask* pf) { return (dynamic_cast<FEBioOpt*>(pf) != 0); }
	FEBioTask* Create(FEModel* pfem) { return new FEBioOpt(pfem); }
};

FEBioOptFactory	febioopt_factory;

//-----------------------------------------------------------------------------
// keep a copy of the FEBioKernel
FEBioKernel*	pFEBio;

//-----------------------------------------------------------------------------
// This function is called from the application and is used to register the plugin classes.
extern "C" DLL_EXPORT void RegisterPlugin(FEBioKernel& febio)
{
	// register the plugin classes
	febio.RegisterClass(&febioopt_factory);

	// store a pointer to the kernel
	pFEBio = &febio;

	// write something to the logfile
	Logfile& log = pFEBio->GetLogfile();
	log.printf("Hello, world!!!\n");
}
