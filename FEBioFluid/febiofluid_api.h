#pragma once

#ifdef WIN32
	#ifdef FECORE_DLL
		#ifdef FEBIOFLUID_EXPORTS
			#define FEBIOFLUID_API __declspec(dllexport)
		#else
			#define FEBIOFLUID_API __declspec(dllimport)
		#endif
	#else
		#define FEBIOFLUID_API
	#endif
#else
	#define FEBIOFLUID_API
#endif
