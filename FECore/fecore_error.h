#pragma once
#include "fecore_api.h"

//-----------------------------------------------------------------------------
// helper functions for reporting errors
FECORE_API bool fecore_error(const char* szerr, ...);
FECORE_API const char* fecore_get_error_string();
