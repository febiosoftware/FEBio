#include "stdafx.h"
#include "version.h"
#include "validate.h"
#include <FECore/log.h>
#include <stdio.h>

unsigned char banner[] = {
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,8,8,8,8,8,8,0,0,0,8,8,8,8,8,8,8,8,0,0,0,8,8,8,8,8,8,8,8,0,0,0,0,8,8,7,0,0,0,0,0,0,8,8,8,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,8,8,8,8,8,8,7,0,0,8,8,8,8,8,8,8,8,7,0,0,8,8,8,8,8,8,8,8,8,7,0,0,0,7,7,0,0,0,0,8,8,8,8,8,8,8,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,7,7,7,7,7,7,7,0,0,8,8,7,7,7,7,7,7,7,0,0,8,8,7,7,7,7,7,8,8,8,7,0,0,0,0,0,0,0,8,8,8,7,7,7,7,8,8,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,8,8,7,0,8,8,0,0,0,0,8,8,7,7,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,8,8,7,0,8,8,7,0,0,8,8,7,7,0,0,0,0,0,0,8,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,8,8,8,7,0,8,8,7,0,0,8,8,7,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,8,8,8,8,8,0,0,0,0,8,8,8,8,8,8,8,0,0,0,0,8,8,8,8,8,8,8,8,8,7,7,0,8,8,7,0,0,8,8,7,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,8,8,8,8,8,7,0,0,0,8,8,8,8,8,8,8,7,0,0,0,8,8,8,8,8,8,8,8,8,7,0,0,8,8,7,0,0,8,8,7,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,7,7,7,7,7,7,0,0,0,8,8,7,7,7,7,7,7,0,0,0,8,8,7,7,7,7,7,8,8,8,0,0,8,8,7,0,0,8,8,7,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,8,8,7,0,8,8,7,0,0,8,8,7,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,8,8,7,0,8,8,7,0,0,0,8,8,0,0,0,0,0,0,8,8,7,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,8,8,8,7,0,8,8,7,0,0,0,8,8,8,0,0,0,0,8,8,8,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,8,8,8,8,8,8,8,8,0,0,0,8,8,8,8,8,8,8,8,8,7,7,0,8,8,7,0,0,0,0,8,8,8,8,8,8,8,8,7,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,8,8,7,0,0,0,0,0,0,0,0,8,8,8,8,8,8,8,8,7,0,0,8,8,8,8,8,8,8,8,7,7,0,0,8,8,7,0,0,0,0,0,7,8,8,8,8,7,7,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,7,0,0,0,0,0,0,0,0,0,7,7,7,7,7,7,7,7,0,0,0,7,7,7,7,7,7,7,7,0,0,0,0,7,7,0,0,0,0,0,0,0,7,7,7,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};


///////////////////////////////////////////////////////////////////////////////
// FUNCTION : Hello
// Prints the FEBio banner to a file
//
/*
#ifdef WIN32

//#include "console.h"
#include "windows.h"

void print_banner()
{
	char sz[] = "                                                                          ";

	Console* pwnd = Console::GetHandle();

	pwnd->Write(sz, 0xF0 );
	pwnd->Draw(banner, 80, 18);

	printf("                                                                           \n");
	printf("      F I N I T E   E L E M E N T S   F O R   B I O M E C H A N I C S      \n");
	printf("                                                                           \n");
	printf("                 --- v e r s i o n - %d . %d . %d", VERSION, SUBVERSION, SUBSUBVERSION);
	if (SVNREVISION) printf(" . %d ---                 \n", SVNREVISION);
	else printf(" ---                 \n");
	printf("                                                                           \n");
	printf("                                                                           \n");
	printf("  Musculoskeletal Research Laboratory                                      \n");
	printf("  University of Utah                                                       \n");
	printf("  http://mrl.sci.utah.edu                                                  \n");
	printf("                                                                           \n");
	printf("  copyright (c) 2006-2013 - All rights reserved                            \n");
	printf("                                                                           \n");
	pwnd->Write(sz, 0xF0 );
	printf("\n\n");
}

#else

void print_banner() {}

#endif
*/
void Hello()
{
	int nlic = GetLicenseKeyStatus();

	clog.printf("===========================================================================\n");
	clog.printf("         ________    _________   _________     __     _________            \n");
	clog.printf("        |        |\\ |        |\\ |        |\\   |  |\\  /         \\\\          \n");
	clog.printf("        |    ____|| |    ____|| |    __  ||   |__|| |    ___    ||         \n");
	clog.printf("        |   |\\___\\| |   |\\___\\| |   |\\_| ||    \\_\\| |   //  \\   ||         \n");
	clog.printf("        |   ||      |   ||      |   || | ||    __   |  ||    |  ||         \n");
	clog.printf("        |   ||__    |   ||__    |   ||_| ||   |  |\\ |  ||    |  ||         \n");
	clog.printf("        |       |\\  |       |\\  |         \\\\  |  || |  ||    |  ||         \n");
	clog.printf("        |    ___||  |    ___||  |    ___   || |  || |  ||    |  ||         \n");
	clog.printf("        |   |\\__\\|  |   |\\__\\|  |   |\\__|  || |  || |  ||    |  ||         \n");
	clog.printf("        |   ||      |   ||      |   ||  |  || |  || |  ||    |  ||         \n");
	clog.printf("        |   ||      |   ||___   |   ||__|  || |  || |   \\\\__/   ||         \n");
	clog.printf("        |   ||      |        |\\ |          || |  || |           ||         \n");
	clog.printf("        |___||      |________|| |__________|| |__||  \\_________//          \n");
	clog.printf("                                                                           \n");
	clog.printf("      F I N I T E   E L E M E N T S   F O R   B I O M E C H A N I C S      \n");
	clog.printf("                                                                           \n");
	clog.printf("                 --- v e r s i o n - %d . %d . %d", VERSION, SUBVERSION, SUBSUBVERSION);
	if (SVNREVISION) clog.printf(" . %d ---                 \n", SVNREVISION);
	else clog.printf(" ---                 \n");
	clog.printf("                                                                           \n");
	clog.printf("                                                                           \n");
	clog.printf("  Musculoskeletal Research Laboratory                                      \n");
	clog.printf("  University of Utah                                                       \n");
	clog.printf("  http://mrl.sci.utah.edu                                                  \n");
	clog.printf("                                                                           \n");
	clog.printf("  copyright (c) 2006-2013 - All rights reserved                            \n");
	if(nlic == 0)
	{
		clog.printf("                                                                              \n");
		clog.printf(" This is the NON-COMMERCIAL version of FEBio or the commercial license        \n");
		clog.printf(" key file could not be found. If you have a key file make sure it is          \n");
		clog.printf(" placed in the same directory as the executable. This version may only        \n");
		clog.printf(" be used for non-commercial purposes as described in the license agreement.   \n");
		clog.printf(" The functionality of this version may be limited compared to the commercial  \n");
		clog.printf(" version. If you wish to obtain a valid commercial license file, please       \n");
		clog.printf(" contact the developers.                                                      \n");
	}
	else if (nlic == 1)
	{
		clog.printf("                                                                           \n");
		clog.printf("  This version is licensed to:                                             \n");
		clog.printf("  \t%s\n", GetLicenseUser());
		clog.printf("  \t%s\n", GetLicenseCompany());
	}
	else if(nlic == 2)
	{
		clog.printf("                                                                           \n");
		clog.printf(" The license file is INVALID. You may continue to use FEBio as an          \n");
		clog.printf(" unlicensed version. This implies that this version of FEBio may only be   \n");
		clog.printf(" used for non-commercial purposes as described in the license agreement and\n");
		clog.printf(" the functionality may be limited. If you wish to obtain a valid license   \n");
		clog.printf(" file or if you think your license is valid, please contact the developers.\n");
	}

	clog.printf("                                                                           \n");
	clog.printf("===========================================================================\n");
	clog.printf("\n\n");
}
