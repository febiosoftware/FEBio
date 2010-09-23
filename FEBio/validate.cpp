#include "stdafx.h"

#ifdef FEBIO_LICENSE

#include <KeyGen\validate.h>
#include <string>
#include <iostream>
using namespace std;

int get_app_path (char *pname, size_t pathsize);

bool validate_err(const char* szerr, FILE* fp = 0)
{
	if (szerr) cerr << szerr << endl;
	if (fp) fclose(fp);
	return false;
}

bool read_str(char* sz, int len, FILE* fp)
{
	if (feof(fp) || ferror(fp)) return false;
	fgets(sz, len, fp);
	if (ferror(fp)) return false;

	char* ch = strrchr(sz, '\n');
	if (ch) *ch = 0;
	ch = strrchr(sz, '\r');
	if (ch) *ch = 0;

	return true;
}

bool validate_license()
{
	// get the application path
	char szpath[4096] = {0};
	get_app_path(szpath, 4095);
	char* ch = strrchr(szpath, '/');
	if (ch == 0) ch = strrchr(szpath, '\\');
	if (ch) *(ch+1) = 0; else szpath[0] = 0;

	// define the filename
	char szfile[4096] = {0};
	sprintf(szfile, "%sfebio.lic", szpath);

	// try to open the license file
	FILE* fp = fopen(szfile, "rt");
	if (fp == 0) return validate_err("Failed opening license file\n\n");

	// read the license file
	char szname[256] = {0};
	if (read_str(szname, 255, fp) == false) return validate_err("Failed reading license file\n\n", fp);

	char szcomp[256] = {0};
	if (read_str(szcomp, 255, fp) == false) return validate_err("Failed reading license file\n\n", fp);

	char szkey[256] = {0};
	if (read_str(szkey, 255, fp) == false) return validate_err("Failed reading license file\n\n", fp);

	// close file
	fclose(fp);

	// validate the key
	validatedType ret = validateKey(szkey, szname, szcomp, "0001");
	if (ret != Passed) return validate_err("Invalid license key\n\n");

	// if we make it here, all is well
	return true;
}

#endif
