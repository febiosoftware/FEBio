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
#include "version.h"
#include "LogStream.h"
#include "febio.h"
#include <stdio.h>

int febio::Hello(LogStream& log)
{
	char szversion[128] = { 0 };
	char* szvernum = getVersionString();
	sprintf(szversion, "  version %s\n", szvernum);
	
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
	log.print("  copyright (c) 2006-2023 - All rights reserved                            \n");
	log.print("                                                                           \n");
	log.print("===========================================================================\n");
	log.print("\n");

	return 0;
}
 