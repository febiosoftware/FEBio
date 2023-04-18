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
#pragma once
#include "febiolib_api.h"
#include <string>
namespace febio {

//! This structures stores the command line options that were input by the user
struct FEBIOLIB_API CMDOPTIONS
{
	enum { MAXFILE = 512 };

	int		ndebug;			//!< debug flag

	bool	bsplash;			//!< show splash screen or not
	bool	bsilent;			//!< run FEBio in silent mode (no output to screen)
	bool	binteractive;		//!< start FEBio interactively

	int		dumpLevel;		//!< requested restart level
	int		dumpStride;		//!< (cold) restart file stride

	char	szfile[MAXFILE];	//!< model input file name
	char	szlog[MAXFILE];	//!< log file name
	char	szplt[MAXFILE];	//!< plot file name
	char	szdmp[MAXFILE];	//!< dump file name
	char	szcnf[MAXFILE];	//!< configuration file
	char	sztask[MAXFILE];	//!< task name
	char	szctrl[MAXFILE];	//!< control file for tasks
	char	szimp[MAXFILE];		//!< import file

	CMDOPTIONS()
	{
		defaults();
	}

	void defaults()
	{
		ndebug = 0;
		bsplash = true;
		bsilent = false;
		binteractive = false;
		dumpLevel = 0;
		dumpStride = 1;

		szfile[0] = 0;
		szlog[0] = 0;
		szplt[0] = 0;
		szdmp[0] = 0;
		szcnf[0] = 0;
		sztask[0] = 0;
		szctrl[0] = 0;
		szimp[0] = 0;
	}
};

	// process string for command line options
	FEBIOLIB_API bool ProcessOptionsString(const std::string& s, CMDOPTIONS& ops);
}
