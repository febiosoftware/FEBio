#include "stdafx.h"
#include "version.h"
#include "validate.h"

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

#ifdef WIN32

#include "console.h"
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
	printf("  copyright (c) 2006-2010 - All rights reserved                            \n");
	printf("                                                                           \n");
	pwnd->Write(sz, 0xF0 );
	printf("\n\n");
}

#else

void print_banner() {}

#endif

void Hello(FILE* fp)
{
	int nlic = GetLicenseKeyStatus();

	fprintf(fp,"===========================================================================\n");
	fprintf(fp,"         ________    _________   _________     __     _________            \n");
	fprintf(fp,"        |        |\\ |        |\\ |        |\\   |  |\\  /         \\\\          \n");
	fprintf(fp,"        |    ____|| |    ____|| |    __  ||   |__|| |    ___    ||         \n");
	fprintf(fp,"        |   |\\___\\| |   |\\___\\| |   |\\_| ||    \\_\\| |   //  \\   ||         \n");
	fprintf(fp,"        |   ||      |   ||      |   || | ||    __   |  ||    |  ||         \n");
	fprintf(fp,"        |   ||__    |   ||__    |   ||_| ||   |  |\\ |  ||    |  ||         \n");
	fprintf(fp,"        |       |\\  |       |\\  |         \\\\  |  || |  ||    |  ||         \n");
	fprintf(fp,"        |    ___||  |    ___||  |    ___   || |  || |  ||    |  ||         \n");
	fprintf(fp,"        |   |\\__\\|  |   |\\__\\|  |   |\\__|  || |  || |  ||    |  ||         \n");
	fprintf(fp,"        |   ||      |   ||      |   ||  |  || |  || |  ||    |  ||         \n");
	fprintf(fp,"        |   ||      |   ||___   |   ||__|  || |  || |   \\\\__/   ||         \n");
	fprintf(fp,"        |   ||      |        |\\ |          || |  || |           ||         \n");
	fprintf(fp,"        |___||      |________|| |__________|| |__||  \\_________//          \n");
	fprintf(fp,"                                                                           \n");
	fprintf(fp,"      F I N I T E   E L E M E N T S   F O R   B I O M E C H A N I C S      \n");
	fprintf(fp,"                                                                           \n");
	fprintf(fp,"                 --- v e r s i o n - %d . %d . %d", VERSION, SUBVERSION, SUBSUBVERSION);
	if (SVNREVISION) fprintf(fp," . %d ---                 \n", SVNREVISION);
	else fprintf(fp," ---                 \n");
	fprintf(fp,"                                                                           \n");
	fprintf(fp,"                                                                           \n");
	fprintf(fp,"  Musculoskeletal Research Laboratory                                      \n");
	fprintf(fp,"  University of Utah                                                       \n");
	fprintf(fp,"  http://mrl.sci.utah.edu                                                  \n");
	fprintf(fp,"                                                                           \n");
	fprintf(fp,"  copyright (c) 2006-2010 - All rights reserved                            \n");
	if(nlic == 0)
	{
		fprintf(fp,"                                                                           \n");
		fprintf(fp," This is an UNLICENSED version of FEBio or the license file could not be   \n");
		fprintf(fp," found. If you have a license file make sure it is placed in the same      \n");
		fprintf(fp," directory as the executable. This version may only be used for            \n");
		fprintf(fp," non-commercial purposes as described in the license agreement. The        \n");
		fprintf(fp," functionality of this version may be limited. If you wish to obtain a     \n");
		fprintf(fp," valid license file, please contact the developers.                        \n");
	}
	else if (nlic == 1)
	{
		fprintf(fp,"                                                                           \n");
		fprintf(fp,"  This version is licensed to:                                              \n");
		fprintf(fp,"  \t%s\n", GetLicenseUser());
		fprintf(fp,"  \t%s\n", GetLicenseCompany());
	}
	else if(nlic == 2)
	{
		fprintf(fp,"                                                                           \n");
		fprintf(fp," The license file is INVALID. You may continue to use FEBio as an          \n");
		fprintf(fp," unlicensed version. This implies that this version of FEBio may only be   \n");
		fprintf(fp," used for non-commercial purposes as described in the license agreement and\n");
		fprintf(fp," the functionality may be limited. If you wish to obtain a valid license   \n");
		fprintf(fp," file or if you think your license is valid, please contact the developers.\n");
	}

	fprintf(fp,"                                                                           \n");
	fprintf(fp,"===========================================================================\n");
	fprintf(fp,"\n\n");
}
