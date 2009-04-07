#include "stdafx.h"

#ifdef WIN32
#include "windows.h"
#endif

#ifdef LINUX
#include <unistd.h>
#endif // LINUX

#ifdef __APPLE__
#include <sys/param.h>
#include <mach-o/dyld.h>
#include <unistd.h>
#endif

//-----------------------------------------------------------------------------
// This function determines *where* the exe actually lives
// This code is from Ian MacArthur from the FLTK forum: fltk.general 
//

int get_app_path (char *pname, size_t pathsize)
{
	long result;

#ifdef LINUX
	/* Oddly, the readlink(2) man page says no NULL is appended. */
	/* So you have to do it yourself, based on the return value: */
	pathsize --; /* Preserve a space to add the trailing NULL */
	result = readlink("/proc/self/exe", pname, pathsize);
	if (result > 0)
	{
		pname[result] = 0; /* add the #@!%ing NULL */

		if ((access(pname, 0) == 0))
		return 0; /* file exists, return OK */
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
		return 0; /* file exists, return OK */
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
	int status = -1;
	char given_path[MAXPATHLEN * 2];
	if (!given_path) return status;

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
	return status;
#endif /* APPLE */

	return -1; /* Path Lookup Failed */
}
