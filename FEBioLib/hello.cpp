#include "stdafx.h"
#include "version.h"
#include "FECore/log.h"
#include "febio.h"
#include <stdio.h>

int febio::Hello(int licenseStatus)
{
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
	felog.printf("  copyright (c) 2006-2018 - All rights reserved                            \n");
	felog.printf("                                                                           \n");

	if (licenseStatus == 0)
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
	else if (licenseStatus == 1)
	{
		// TODO: print something useful here. E.g. when the license will expire or something?
	}
	else if (licenseStatus < 0)
	{
		felog.printf("                                                                           \n");
		felog.printf(" The license file is INVALID. You may continue to use FEBio as an          \n");
		felog.printf(" unlicensed version. This implies that this version of FEBio may only be   \n");
		felog.printf(" used for non-commercial purposes as described in the license agreement and\n");
		felog.printf(" the functionality may be limited. If you wish to obtain a valid license   \n");
		felog.printf(" file or if you think your license is valid, please contact the developers.\n");
	}
	else if (licenseStatus == 2)
	{
		felog.printf("                                                                           \n");
		felog.printf(" This is a TRIAL license.                                                  \n");
	}

	felog.printf("                                                                           \n");
	felog.printf("===========================================================================\n");
	felog.printf("\n\n");

	return 0;
}
