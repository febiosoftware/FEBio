#pragma once

#ifdef WIN32
	#ifdef FECORE_DLL
		#ifdef FECORE_EXPORTS
			#define FECORE_API __declspec(dllexport)
		#else
			#define FECORE_API __declspec(dllimport)
		#endif
	#else
		#define FECORE_API
	#endif
#else
	#define FECORE_API
#endif
