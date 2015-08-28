
#ifndef KEYGEN_COMMON_H
#define KEYGEN_COMMON_H 1

//
// common.h
//
// Copyright (c) 2010 - University of Utah Software Development Center
//

#define KEY_SIZE 32   // KEY_SIZE should not be under 14 otherwise the generate/validate functions will need to be changed
#define MOD_SALT 13   // MOD_SALT should not be over 26 otherwise the generate algorithm will need to be changed

#ifdef _WIN32 // _WIN32 is defined by all Windows 32 compilers, but not by others.
#  define EXPORT __declspec(dllexport)
#else
#  define EXPORT
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Not sure if I like 'ZERO_CHAR' or if we should just put '0' in the code... 
// Certainly we should not put the number 48.
//
// 48 - 57 = 0 - 9
// 65 - 90 = A - Z
const char ZERO_CHAR = '0';
const char NINE_CHAR = '9';
const char A_CHAR    = 'A';
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <string>

std::string hashNameAndCompany( const std::string & name, const std::string & company );

#endif // KEYGEN_COMMON_H
