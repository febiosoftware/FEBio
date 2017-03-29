#pragma once

#ifdef WIN32
	#ifdef FECORE_DLL
		#ifdef FECORE_API
			#define FECOREDLL_EXPORT __declspec(dllexport)
		#else
			#define FECOREDLL_EXPORT __declspec(dllimport)
		#endif
	#else
		#define FECOREDLL_EXPORT
	#endif
#else
	#define FECOREDLL_EXPORT
#endif
