#pragma once

#ifdef WIN32
#include <float.h>
#define ISNAN(x) _isnan(x)
#endif

#ifdef LINUX
	#ifdef CENTOS
		#define ISNAN(x) isnan(x)
	#else
		#define ISNAN(x) std::isnan(x)
	#endif
#endif

#ifdef __APPLE__
#include <math.h>
#define ISNAN(x) isnan(x)
#endif

#ifdef WIN32
extern "C" int __cdecl omp_get_num_threads(void);
extern "C" int __cdecl omp_get_thread_num(void);
#else
extern "C" int omp_get_num_threads(void);
extern "C" int omp_get_thread_num(void);
#endif
