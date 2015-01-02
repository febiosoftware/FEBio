#include "stdafx.h"
#include "version.h"
#include "validate.h"
#include "FECore/log.h"
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
	printf("  copyright (c) 2006-2015 - All rights reserved                            \n");
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

	felog.printf("===========================================================================\n");
	felog.printf("         ________    _________   _________     __     _________            \n");
	felog.printf("        |        |\\ |        |\\ |        |\\   |  |\\  /         \\\\          \n");
	felog.printf("        |    ____|| |    ____|| |    __  ||   |__|| |    ___    ||         \n");
	felog.printf("        |   |\\___\\| |   |\\___\\| |   |\\_| ||    \\_\\| |   //  \\   ||         \n");
	felog.printf("        |   ||      |   ||      |   || | ||    __   |  ||    |  ||         \n");
	felog.printf("        |   ||__    |   ||__    |   ||_| ||   |  |\\ |  ||    |  ||         \n");
	felog.printf("        |       |\\  |       |\\  |         \\\\  |  || |  ||    |  ||         \n");
	felog.printf("        |    ___||  |    ___||  |    ___   || |  || |  ||    |  ||         \n");
	felog.printf("        |   |\\__\\|  |   |\\__\\|  |   |\\__|  || |  || |  ||    |  ||         \n");
	felog.printf("        |   ||      |   ||      |   ||  |  || |  || |  ||    |  ||         \n");
	felog.printf("        |   ||      |   ||___   |   ||__|  || |  || |   \\\\__/   ||         \n");
	felog.printf("        |   ||      |        |\\ |          || |  || |           ||         \n");
	felog.printf("        |___||      |________|| |__________|| |__||  \\_________//          \n");
	felog.printf("                                                                           \n");
	felog.printf("      F I N I T E   E L E M E N T S   F O R   B I O M E C H A N I C S      \n");
	felog.printf("                                                                           \n");
	felog.printf("                 --- v e r s i o n - %d . %d . %d", VERSION, SUBVERSION, SUBSUBVERSION);
	if (SVNREVISION) felog.printf(" . %d ---                 \n", SVNREVISION);
	else felog.printf(" ---                 \n");
	felog.printf("                                                                           \n");
	felog.printf("                                                                           \n");
	felog.printf("  Musculoskeletal Research Laboratory                                      \n");
	felog.printf("  University of Utah                                                       \n");
	felog.printf("  http://febio.org                                                         \n");
	felog.printf("                                                                           \n");
	felog.printf("  FEBio is a registered trademark.                                         \n");
	felog.printf("  copyright (c) 2006-2015 - All rights reserved                            \n");
	if(nlic == 0)
	{
		felog.printf("                                                                              \n");
		felog.printf(" This is the NON-COMMERCIAL version of FEBio or the commercial license        \n");
		felog.printf(" key file could not be found. If you have a key file make sure it is          \n");
		felog.printf(" placed in the same directory as the executable. This version may only        \n");
		felog.printf(" be used for non-commercial purposes as described in the license agreement.   \n");
		felog.printf(" The functionality of this version may be limited compared to the commercial  \n");
		felog.printf(" version. If you wish to obtain a valid commercial license file, please       \n");
		felog.printf(" contact the developers.                                                      \n");
	}
	else if (nlic == 1)
	{
		felog.printf("                                                                           \n");
		felog.printf("  This version is licensed to:                                             \n");
		felog.printf("  \t%s\n", GetLicenseUser());
		felog.printf("  \t%s\n", GetLicenseCompany());
	}
	else if(nlic == 2)
	{
		felog.printf("                                                                           \n");
		felog.printf(" The license file is INVALID. You may continue to use FEBio as an          \n");
		felog.printf(" unlicensed version. This implies that this version of FEBio may only be   \n");
		felog.printf(" used for non-commercial purposes as described in the license agreement and\n");
		felog.printf(" the functionality may be limited. If you wish to obtain a valid license   \n");
		felog.printf(" file or if you think your license is valid, please contact the developers.\n");
	}

	felog.printf("                                                                           \n");
	felog.printf("===========================================================================\n");
	felog.printf("\n\n");
}
