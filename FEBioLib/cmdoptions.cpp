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
#include "cmdoptions.h"
#include "febio.h"
#include <stdlib.h>

std::vector< std::string > split_string(const std::string& s)
{
	std::vector< std::string > args;
	std::string t;
	bool inquote = false;
	for (int n = 0; n < s.size(); ++n)
	{
		char c = s[n];
		if (isspace(c) && (inquote== false))
		{
			if (t.empty() == false)
			{
				args.push_back(t);
				t.clear();
			}
		}
		else if (c == '\"')
		{
			if (inquote == false) inquote = true;
			else inquote = false;
		}
		else t.push_back(c);
	}
	if (t.empty() == false)
	{
		args.push_back(t);
	}
	return args;
}

// Copy a string into a fixed-size buffer. Reports an error if the string does not fit,
// instead of overflowing the buffer.
static bool copyString(char* dst, size_t size, const std::string& src)
{
	if (src.size() >= size)
	{
		fprintf(stderr, "FATAL ERROR: file name is too long (%d characters; maximum is %d):\n%s\n", (int)src.size(), (int)size - 1, src.c_str());
		dst[0] = 0;
		return false;
	}
	memcpy(dst, src.c_str(), src.size() + 1);
	return true;
}

// Remove the extension from a file name. Only a '.' in the file name itself counts,
// not one in a directory name (e.g. "/path/2. Model/file.feb").
static std::string removeExtension(const std::string& file)
{
	size_t dot = file.rfind('.');
	size_t slash = file.find_last_of("/\\");
	if ((dot != std::string::npos) && ((slash == std::string::npos) || (dot > slash)))
		return file.substr(0, dot);
	return file;
}

bool febio::ProcessOptionsString(const std::string& s, CMDOPTIONS& ops)
{
	// split the string
	std::vector< std::string > args = split_string(s);

	// set default options
	ops.ndebug = 0;
	ops.bsplash = true;
	ops.bsilent = false;
	ops.binteractive = true;

	// these flags indicate whether the corresponding file name
	// was defined on the command line. Otherwise, a default name will be generated.
	bool blog = false;
	bool bplt = false;
	bool bdmp = false;
	bool brun = true;

	// initialize file names
	ops.szfile[0] = 0;
	ops.szplt[0] = 0;
	ops.szlog[0] = 0;
	ops.szdmp[0] = 0;
	ops.sztask[0] = 0;
	ops.szctrl[0] = 0;
	ops.szimp[0] = 0;

	// set initial configuration file name
	if (ops.szcnf[0] == 0)
	{
		char szpath[1024] = { 0 };
		febio::get_app_path(szpath, 1023);
		snprintf(ops.szcnf, sizeof(ops.szcnf), "%sfebio.xml", szpath);
	}

	// loop over the arguments
	int nargs = args.size();
	for (int i = 1; i < nargs; ++i)
	{
		const char* sz = args[i].c_str();

		if (strcmp(sz, "-r") == 0)
		{
			if (ops.sztask[0] != 0) { fprintf(stderr, "-r is incompatible with other command line option.\n"); return false; }
			copyString(ops.sztask, sizeof(ops.sztask), "restart");
			copyString(ops.szctrl, sizeof(ops.szctrl), args[++i].c_str());
			ops.binteractive = false;
		}
		else if (strcmp(sz, "-noappend") == 0)
		{
			ops.bappendFiles = false;
		}
		else if (strcmp(sz, "-d") == 0)
		{
			if (ops.sztask[0] != 0) { fprintf(stderr, "-d is incompatible with other command line option.\n"); return false; }
			copyString(ops.sztask, sizeof(ops.sztask), "diagnose");
			copyString(ops.szctrl, sizeof(ops.szctrl), args[++i].c_str());
			ops.binteractive = false;
		}
		else if (strcmp(sz, "-p") == 0)
		{
			bplt = true;
			copyString(ops.szplt, sizeof(ops.szplt), args[++i].c_str());
		}
		else if (strncmp(sz, "-dump", 5) == 0)
		{
			ops.dumpLevel = FE_DUMP_MAJOR_ITRS;
			if (sz[5] == '=') ops.dumpLevel = atoi(sz + 6);
			if ((ops.dumpLevel < 0) || (ops.dumpLevel > 3))
			{
				fprintf(stderr, "FATAL ERROR: invalid restart level.\n");
				return false;
			}

			if (i < nargs - 1)
			{
				const char* szi = args[i + 1].c_str();
				if (szi[0] != '-')
				{
					// assume this is the name of the dump file
					copyString(ops.szdmp, sizeof(ops.szdmp), args[++i].c_str());
					bdmp = true;
				}
			}
		}
		else if (strcmp(sz, "-o") == 0)
		{
			blog = true;
			copyString(ops.szlog, sizeof(ops.szlog), args[++i].c_str());
		}
		else if (strcmp(sz, "-i") == 0)
		{
			if (i >= nargs - 1) { fprintf(stderr, "FATAL ERROR: missing file name after -i\n"); return false; }
			sz = args[++i].c_str();
			const char* szext = strrchr(sz, '.');
			if (szext == 0)
			{
				// we assume a default extension of .feb if none is provided
				snprintf(ops.szfile, sizeof(ops.szfile), "%s.feb", sz);
			}
			else if (!copyString(ops.szfile, sizeof(ops.szfile), sz)) return false;
			ops.binteractive = false;
		}
		else if (strcmp(sz, "-s") == 0)
		{
			if (ops.sztask[0] != 0) { fprintf(stderr, "-s is incompatible with other command line option.\n"); return false; }
			copyString(ops.sztask, sizeof(ops.sztask), "optimize");
			copyString(ops.szctrl, sizeof(ops.szctrl), args[++i].c_str());
		}
		else if ((strcmp(sz, "-g") == 0) || (strcmp(sz, "-g1") == 0))
		{
			ops.ndebug = 1;
		}
		else if (strcmp(sz, "-g2") == 0)
		{
			ops.ndebug = 2;
		}
		else if (strcmp(sz, "-nosplash") == 0)
		{
			// don't show the welcome message
			ops.bsplash = false;
		}
		else if (strcmp(sz, "-silent") == 0)
		{
			// no output to screen
			ops.bsilent = true;
		}
		else if (strcmp(sz, "-cnf") == 0)	// obsolete: use -config instead
		{
			copyString(ops.szcnf, sizeof(ops.szcnf), args[++i].c_str());
		}
		else if (strcmp(sz, "-config") == 0)
		{
			copyString(ops.szcnf, sizeof(ops.szcnf), args[++i].c_str());
		}
		else if (strcmp(sz, "-noconfig") == 0)
		{
			ops.szcnf[0] = 0;
		}
		else if (strncmp(sz, "-task", 5) == 0)
		{
			if (sz[5] != '=') { fprintf(stderr, "command line error when parsing task\n"); return false; }
			copyString(ops.sztask, sizeof(ops.sztask), sz + 6);

			if (i < nargs - 1)
			{
				const char* szi = args[i + 1].c_str();
				if (szi[0] != '-')
				{
					// assume this is a control file for the specified task
					copyString(ops.szctrl, sizeof(ops.szctrl), args[++i].c_str());
				}
			}
		}
		else if (strcmp(sz, "-import") == 0)
		{
			copyString(ops.szimp, sizeof(ops.szimp), args[++i].c_str());
		}
		else if (sz[0] == '-')
		{
			fprintf(stderr, "FATAL ERROR: Invalid command line option '%s'.\n", sz);
			return false;
		}
		else
		{
			// if no input file is given yet, we'll assume this is the input file
			if (ops.szfile[0] == 0)
			{
				const char* szext = strrchr(sz, '.');
				if (szext == 0)
				{
					// we assume a default extension of .feb if none is provided
					snprintf(ops.szfile, sizeof(ops.szfile), "%s.feb", sz);
				}
				else
				{
					if (!copyString(ops.szfile, sizeof(ops.szfile), sz)) return false;
				}
				ops.binteractive = false;
			}
			else
			{
				fprintf(stderr, "FATAL ERROR: Invalid command line option '%s'\n", sz);
				return false;
			}
		}
	}

	// do some sanity checks
	if (strcmp(ops.sztask, "optimize") == 0)
	{
		// make sure we have an input file
		if (ops.szfile[0] == 0)
		{
			fprintf(stderr, "FATAL ERROR: no model input file was defined (use -i to define the model input file)\n\n");
			return false;
		}
	}

	// if no task is defined, we assume a std solve is wanted
	if (ops.sztask[0] == 0) copyString(ops.sztask, sizeof(ops.sztask), "solve");

	// derive the other filenames
	// NOTE: std::string is used here so that long file paths are not truncated
	// (or overflow a fixed-size buffer).
	if (ops.szfile[0] || ops.szctrl[0])
	{
		std::string sbase = removeExtension(ops.szfile[0] ? ops.szfile : ops.szctrl);

		std::string slogbase = sbase;
		if (ops.szfile[0] && ops.szctrl[0]) slogbase = removeExtension(ops.szctrl);

		if (!blog && !copyString(ops.szlog, sizeof(ops.szlog), slogbase + ".log")) return false;
		if (!bplt && !copyString(ops.szplt, sizeof(ops.szplt), sbase + ".xplt")) return false;
		if (!bdmp && !copyString(ops.szdmp, sizeof(ops.szdmp), sbase + ".dmp")) return false;
	}

	return true;
}
