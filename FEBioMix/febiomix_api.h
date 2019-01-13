#pragma once

#ifdef WIN32
	#ifdef FECORE_DLL
		#ifdef FEBIOMIX_EXPORTS
			#define FEBIOMIX_API __declspec(dllexport)
		#else
			#define FEBIOMIX_API __declspec(dllimport)
		#endif
	#else
		#define FEBIOMIX_API
	#endif
#else
	#define FEBIOMIX_API
#endif
