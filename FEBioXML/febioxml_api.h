#pragma once

#ifdef WIN32
	#ifdef FECORE_DLL
		#ifdef FEBIOXML_EXPORTS
			#define FEBIOXML_API __declspec(dllexport)
		#else
			#define FEBIOXML_API __declspec(dllimport)
		#endif
	#else
		#define FEBIOXML_API
	#endif
#else
	#define FEBIOXML_API
#endif
