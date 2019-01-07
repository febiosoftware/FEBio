#pragma once

#ifdef WIN32
	#ifdef USE_DLL
		#ifdef FEBIOLIB_EXPORTS
			#define FEBIOLIB_API __declspec(dllexport)
		#else
			#define FEBIOLIB_API __declspec(dllimport)
		#endif
	#else
		#define FEBIOLIB_API
	#endif
#else
	#define FEBIOLIB_API
#endif
