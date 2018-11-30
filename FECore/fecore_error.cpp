#include "stdafx.h"
#include "fecore_error.h"
#include "FECoreKernel.h"
#include <stdarg.h>

//-----------------------------------------------------------------------------
//! Helper function for reporting errors
bool fecore_error(const char* sz, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// make the message
	char szerr[512] = { 0 };
	va_start(args, sz);
	vsprintf(szerr, sz, args);
	va_end(args);

	// TODO: Perhaps I should report it to the logfile?
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	fecore.SetErrorString(szerr);

	return false;
}

//-----------------------------------------------------------------------------
const char* fecore_get_error_string()
{
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	return fecore.GetErrorString();
}
