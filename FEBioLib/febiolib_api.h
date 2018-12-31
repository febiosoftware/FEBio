#pragma once

#ifdef WIN32
	#ifdef FEBIOLIB_DLL
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
