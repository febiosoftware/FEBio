#include "stdafx.h"
#include "version.h"
#include <FECore/log.h>
#include "febio.h"
#include <stdio.h>

int febio::Hello()
{
	felog.printf("===========================================================================\n");
	felog.printf("         ________    _________   _______       __     _________            \n");
	felog.printf("        |        |\\ |        |\\ |       \\\\    |  |\\  /         \\\\          \n");
	felog.printf("        |    ____|| |    ____|| |    __  ||   |__|| |    ___    ||         \n");
	felog.printf("        |   |\\___\\| |   |\\___\\| |   |\\_| ||    \\_\\| |   //  \\   ||         \n");
	felog.printf("        |   ||__    |   ||__    |   ||_| ||   |  |\\ |  ||    |  ||         \n");
	felog.printf("        |       |\\  |       |\\  |         \\\\  |  || |  ||    |  ||         \n");
	felog.printf("        |    ___||  |    ___||  |    ___   || |  || |  ||    |  ||         \n");
	felog.printf("        |   |\\__\\|  |   |\\__\\|  |   |\\__|  || |  || |  ||    |  ||         \n");
	felog.printf("        |   ||      |   ||___   |   ||__|  || |  || |   \\\\__/   ||         \n");
	felog.printf("        |   ||      |        |\\ |          || |  || |           ||         \n");
	felog.printf("        |___||      |________|| |_________//  |__||  \\_________//          \n");
	felog.printf("                                                                           \n");
	felog.printf("      F I N I T E   E L E M E N T S   F O R   B I O M E C H A N I C S      \n");
	felog.printf("                                                                           \n");
	felog.printf("  version %d.%d.%d                                                         \n", VERSION, SUBVERSION, SUBSUBVERSION);
	felog.printf("  FEBio is a registered trademark.                                         \n");
	felog.printf("  copyright (c) 2006-2019 - All rights reserved                            \n");
	felog.printf("                                                                           \n");
	felog.printf("===========================================================================\n");
	felog.printf("\n");

	return 0;
}
 