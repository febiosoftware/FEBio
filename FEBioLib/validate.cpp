#include "stdafx.h"
#include "febio.h"
#include "validate.h"
#include "febio.h"
#include <string.h>

#ifndef FEBIOLM
int GetLicenseKeyStatus(const char* licenseKey)
{
	return 0;
}
#else
//#include <febiolm.h>

#endif

std::string LoadLicenseKey()
{
	// Get the location of the executable
	char szpath[1024] = { 0 };
	febio::get_app_path(szpath, 1023);

	// append the file name
	strcat(szpath, "license.txt");

	// try to open the file
	std::string licenseKey;
	FILE* fp = fopen(szpath, "rt");
	if (fp != NULL)
	{
		char szbuf[128] = {0};
		fgets(szbuf, 127, fp);
		licenseKey = szbuf;
		fclose(fp);
	}
	return licenseKey;
}
