/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "febio.h"

#ifdef WIN32
#include "windows.h"
#endif

#ifdef LINUX
#include <unistd.h>
#include <string.h>
#endif // LINUX

#ifdef __APPLE__
#include <sys/param.h>
#include <mach-o/dyld.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#endif

//-----------------------------------------------------------------------------
// This function determines *where* the exe actually lives
// This code is from Ian MacArthur from the FLTK forum: fltk.general 
//

int febio::get_app_path (char *pname, size_t pathsize)
{
	long result;
	int status = -1;
#ifdef LINUX
	/* Oddly, the readlink(2) man page says no NULL is appended. */
	/* So you have to do it yourself, based on the return value: */
	pathsize --; /* Preserve a space to add the trailing NULL */
	result = readlink("/proc/self/exe", pname, pathsize);
	if (result > 0)
	{
		pname[result] = 0; /* add the #@!%ing NULL */

		if ((access(pname, 0) == 0))
		status = 0; /* file exists, return OK */
		          /*else name doesn't seem to exist, return FAIL (falls through) */
	}
#endif /* LINUX */

#ifdef WIN32
	result = GetModuleFileNameA(GetModuleHandle(NULL), (LPCH) pname, pathsize);
	if (result > 0)
	{
		/* fix up the dir slashes... */
		int len = strlen(pname);
		int idx;
		for (idx = 0; idx < len; idx++)
		{
			if (pname[idx] == '\\') pname[idx] = '/';
		}
		status = 0; /* file exists, return OK */
		          /*else name doesn't seem to exist, return FAIL (falls through) */
	}
#endif /* WIN32 */


#ifdef __APPLE__ /* assume this is OSX */
/*
        from http://www.hmug.org/man/3/NSModule.html

        extern int _NSGetExecutablePath(char *buf, unsigned long *bufsize);

        _NSGetExecutablePath  copies  the  path  of the executable
        into the buffer and returns 0 if the path was successfully
        copied  in the provided buffer. If the buffer is not large
        enough, -1 is returned and the  expected  buffer  size  is
        copied  in  *bufsize.  Note that _NSGetExecutablePath will
        return "a path" to the executable not a "real path" to the
        executable.  That  is  the path may be a symbolic link and
        not the real file. And with  deep  directories  the  total
        bufsize needed could be more than MAXPATHLEN.
*/
	char given_path[MAXPATHLEN * 2];

	pathsize = MAXPATHLEN * 2;
	uint32_t nsize = (uint32_t) (pathsize);
	result = _NSGetExecutablePath(given_path, &nsize);
	if (result == 0)
	{ /* OK, we got something - now try and resolve the real path...*/
		if (realpath(given_path, pname) != NULL)
		{
			if ((access(pname, 0) == 0))
			status = 0; /* file exists, return OK */
		}
	}
#endif /* APPLE */

	if (status == 0)
	{
		// remove the application's name
		char* ch = strrchr(pname, '\\');
		if (ch == 0) ch = strrchr(pname, '/');
		if (ch) ch[1] = 0;
	}

	return status;
}
