#pragma once

#ifdef WIN32
	#ifdef FECORE_DLL
		#ifdef FEBIOMECH_EXPORTS
			#define FEBIOMECH_API __declspec(dllexport)
		#else
			#define FEBIOMECH_API __declspec(dllimport)
		#endif
	#else
		#define FEBIOMECH_API
	#endif
#else
	#define FEBIOMECH_API
#endif
