#include "stdafx.h"
#include "version.h"
#include <FECore/Logfile.h>
#include "febio.h"
#include <stdio.h>

int febio::Hello(LogStream& log)
{
	char szversion[128] = { 0 };
	sprintf(szversion, "  version %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
	log.print("===========================================================================\n");
	log.print("         ________    _________   _______       __     _________            \n");
	log.print("        |        |\\ |        |\\ |       \\\\    |  |\\  /         \\\\          \n");
	log.print("        |    ____|| |    ____|| |    __  ||   |__|| |    ___    ||         \n");
	log.print("        |   |\\___\\| |   |\\___\\| |   |\\_| ||    \\_\\| |   //  \\   ||         \n");
	log.print("        |   ||__    |   ||__    |   ||_| ||   |  |\\ |  ||    |  ||         \n");
	log.print("        |       |\\  |       |\\  |         \\\\  |  || |  ||    |  ||         \n");
	log.print("        |    ___||  |    ___||  |    ___   || |  || |  ||    |  ||         \n");
	log.print("        |   |\\__\\|  |   |\\__\\|  |   |\\__|  || |  || |  ||    |  ||         \n");
	log.print("        |   ||      |   ||___   |   ||__|  || |  || |   \\\\__/   ||         \n");
	log.print("        |   ||      |        |\\ |          || |  || |           ||         \n");
	log.print("        |___||      |________|| |_________//  |__||  \\_________//          \n");
	log.print("                                                                           \n");
	log.print("      F I N I T E   E L E M E N T S   F O R   B I O M E C H A N I C S      \n");
	log.print("                                                                           \n");
	log.print(szversion);
	log.print("  FEBio is a registered trademark.                                         \n");
	log.print("  copyright (c) 2006-2019 - All rights reserved                            \n");
	log.print("                                                                           \n");
	log.print("===========================================================================\n");
	log.print("\n");

	return 0;
}
 