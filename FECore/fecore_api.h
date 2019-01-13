#pragma once

#ifdef WIN32
	#ifdef USE_DLL
		#define FECORE_API __declspec(dllexport)
	#else
		#define FECORE_API
	#endif
#else
	#define FECORE_API
#endif
