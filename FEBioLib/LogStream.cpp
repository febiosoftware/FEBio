#include "stdafx.h"
#include "LogStream.h"
#include <stdarg.h>
#include <stdio.h>

void LogStream::printf(const char* sz, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// make the message
	char sztxt[1024] = { 0 };
	va_start(args, sz);
	vsprintf(sztxt, sz, args);
	va_end(args);

	print(sztxt);
}
