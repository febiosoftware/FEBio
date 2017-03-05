#pragma once

#ifdef WIN32
	#ifdef FECORE_DLL
		#ifdef FECORE_API
			#define FECORE_EXPORT __declspec(dllexport)
		#else
			#define FECORE_EXPORT __declspec(dllimport)
		#endif
	#else
		#define FECORE_EXPORT
	#endif
#else
	#define FECORE_EXPORT
#endif
