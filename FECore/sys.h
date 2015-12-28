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
