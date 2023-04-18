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
		sprintf(ops.szcnf, "%sfebio.xml", szpath);
	}

	// loop over the arguments
	int nargs = args.size();
	for (int i = 1; i < nargs; ++i)
	{
		const char* sz = args[i].c_str();

		if (strcmp(sz, "-r") == 0)
		{
			if (ops.sztask[0] != 0) { fprintf(stderr, "-r is incompatible with other command line option.\n"); return false; }
			strcpy(ops.sztask, "restart");
			strcpy(ops.szctrl, args[++i].c_str());
			ops.binteractive = false;
		}
		else if (strcmp(sz, "-d") == 0)
		{
			if (ops.sztask[0] != 0) { fprintf(stderr, "-d is incompatible with other command line option.\n"); return false; }
			strcpy(ops.sztask, "diagnose");
			strcpy(ops.szctrl, args[++i].c_str());
			ops.binteractive = false;
		}
		else if (strcmp(sz, "-p") == 0)
		{
			bplt = true;
			strcpy(ops.szplt, args[++i].c_str());
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
					strcpy(ops.szdmp, args[++i].c_str());
					bdmp = true;
				}
			}
		}
		else if (strcmp(sz, "-o") == 0)
		{
			blog = true;
			strcpy(ops.szlog, args[++i].c_str());
		}
		else if (strcmp(sz, "-i") == 0)
		{
			++i;
			const char* szext = strrchr(sz, '.');
			if (szext == 0)
			{
				// we assume a default extension of .feb if none is provided
				sprintf(ops.szfile, "%s.feb", sz);
			}
			else strcpy(ops.szfile, sz);
			ops.binteractive = false;
		}
		else if (strcmp(sz, "-s") == 0)
		{
			if (ops.sztask[0] != 0) { fprintf(stderr, "-s is incompatible with other command line option.\n"); return false; }
			strcpy(ops.sztask, "optimize");
			strcpy(ops.szctrl, args[++i].c_str());
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
			strcpy(ops.szcnf, args[++i].c_str());
		}
		else if (strcmp(sz, "-config") == 0)
		{
			strcpy(ops.szcnf, args[++i].c_str());
		}
		else if (strcmp(sz, "-noconfig") == 0)
		{
			ops.szcnf[0] = 0;
		}
		else if (strncmp(sz, "-task", 5) == 0)
		{
			if (sz[5] != '=') { fprintf(stderr, "command line error when parsing task\n"); return false; }
			strcpy(ops.sztask, sz + 6);

			if (i < nargs - 1)
			{
				const char* szi = args[i + 1].c_str();
				if (szi[0] != '-')
				{
					// assume this is a control file for the specified task
					strcpy(ops.szctrl, args[++i].c_str());
				}
			}
		}
		else if (strcmp(sz, "-import") == 0)
		{
			strcpy(ops.szimp, args[++i].c_str());
		}
		else if (sz[0] == '-')
		{
			fprintf(stderr, "FATAL ERROR: Invalid command line option.\n");
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
					sprintf(ops.szfile, "%s.feb", sz);
				}
				else
				{
					strcpy(ops.szfile, sz);
				}
				ops.binteractive = false;
			}
			else
			{
				fprintf(stderr, "FATAL ERROR: Invalid command line option\n");
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
	if (ops.sztask[0] == 0) strcpy(ops.sztask, "solve");

	// derive the other filenames
	if (ops.szfile[0])
	{
		char szbase[256]; strcpy(szbase, ops.szfile);
		char* ch = strrchr(szbase, '.');
		if (ch) *ch = 0;

		char szlogbase[256];
		if (ops.szctrl[0])
		{
			strcpy(szlogbase, ops.szctrl);
			ch = strrchr(szlogbase, '.');
			if (ch) *ch = 0;
		}
		else strcpy(szlogbase, szbase);

		if (!blog) sprintf(ops.szlog, "%s.log", szlogbase);
		if (!bplt) sprintf(ops.szplt, "%s.xplt", szbase);
		if (!bdmp) sprintf(ops.szdmp, "%s.dmp", szbase);
	}
	else if (ops.szctrl[0])
	{
		char szbase[256]; strcpy(szbase, ops.szfile);
		strcpy(szbase, ops.szctrl);
		char* ch = strrchr(szbase, '.');
		if (ch) *ch = 0;

		if (!blog) sprintf(ops.szlog, "%s.log", szbase);
		if (!bplt) sprintf(ops.szplt, "%s.xplt", szbase);
		if (!bdmp) sprintf(ops.szdmp, "%s.dmp", szbase);
	}

	return true;
}
