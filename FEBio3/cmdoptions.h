#pragma once

//! This structures stores the command line options that were input by the user
struct CMDOPTIONS
{
	enum { MAXFILE = 512 };

	bool		bdebug;			//!< debug flag

	bool	bsplash;			//!< show splash screen or not
	bool	bsilent;			//!< run FEBio in silent mode (no output to screen)
	bool	binteractive;		//!< start FEBio interactively

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
		bdebug = false;
		bsplash = true;
		bsilent = false;
		binteractive = false;

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
