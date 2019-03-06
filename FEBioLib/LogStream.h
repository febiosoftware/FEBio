#pragma once
#include "febiolib_api.h"

//-----------------------------------------------------------------------------
// class used to create an abstract interface to a screen
class FEBIOLIB_API LogStream
{
public:
	LogStream() {}
	virtual ~LogStream() {}

	// override function to print
	virtual void print(const char* sz) = 0;

	// function to print variable output
	void printf(const char* sz, ...);

	// flush the stream
	virtual void flush() {}
};
