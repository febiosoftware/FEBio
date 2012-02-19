#include "stdafx.h"

#include "validate.h"
#include <string>
#include <iostream>
#include <cstring>
using namespace std;

// load the key generator
#ifdef FEBIO_LICENSE
	#include "KeyGen/validate.h"
#endif

//-----------------------------------------------------------------------------
// FEBio license-key manager
//
class FEBioLicenseKey
{
	// license key status
	enum {
		FEBIO_NO_LICENSE,
		FEBIO_VALID_LICENSE,
		FEBIO_INVALID_LICENSE
	};

public:
	// return the one-and-only key manager
	static FEBioLicenseKey* GetInstance();

	// load a key from file
	void Load();

	// get the license status
	int Status() { return m_nstatus; }

	// return the user
	const char* User();

	// return the company
	const char* Company();

protected:
	void validate_error(FILE* fp = 0);
	
private:
	// hide constructors
	FEBioLicenseKey() { m_nstatus = FEBIO_NO_LICENSE; }
	FEBioLicenseKey(const FEBioLicenseKey& key) {}

	// the one-and-only key
	static FEBioLicenseKey*	m_pKey;

	int		m_nstatus;
	string	m_key;
	string	m_user;
	string	m_comp;
};

//-----------------------------------------------------------------------------
// The one-and-only license key
FEBioLicenseKey* FEBioLicenseKey::m_pKey = 0;

//-----------------------------------------------------------------------------
FEBioLicenseKey* FEBioLicenseKey::GetInstance()
{
	if (m_pKey == 0) m_pKey = new FEBioLicenseKey;
	return m_pKey;
}

//-----------------------------------------------------------------------------
// helper function to read a string from a file
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

//-----------------------------------------------------------------------------
// function to return the application path
extern int get_app_path (char *pname, size_t pathsize);

//-----------------------------------------------------------------------------
// Try to load a license key
void FEBioLicenseKey::Load()
{
#ifdef FEBIO_LICENSE
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
	if (fp == 0) { m_nstatus = FEBIO_NO_LICENSE; return; }

	// read the license key
	char szkey[256] = {0};
	if (read_str(szkey, 255, fp) == false) { validate_error(fp); return; }
	m_key = szkey;

	// read the user name
	char szname[256] = {0};
	if (read_str(szname, 255, fp) == false) { validate_error(fp); return; }
	m_user = szname;

	// read the company name
	char szcomp[256] = {0};
	if (read_str(szcomp, 255, fp) == false) { validate_error(fp); return; }
	m_comp = szcomp;

	// close file
	fclose(fp);

	// validate the key
	ValidationResult ret = validateKey(szkey, szname, szcomp, "0001");
	if (ret != Passed) return validate_error();

	// if we make it here, all is well
	m_nstatus = FEBIO_VALID_LICENSE;
#endif
}

//-----------------------------------------------------------------------------
// helper function to process an error while reading the license file
void FEBioLicenseKey::validate_error(FILE* fp)
{
	if (fp) fclose(fp);
	m_nstatus = FEBIO_INVALID_LICENSE;
	m_key.clear();
	m_user.clear();
	m_comp.clear();
}

//-----------------------------------------------------------------------------
const char* FEBioLicenseKey::User()
{
	if (m_nstatus != FEBIO_VALID_LICENSE) return 0;
	return m_user.c_str();
}

//-----------------------------------------------------------------------------
const char* FEBioLicenseKey::Company()
{
	if (m_nstatus != FEBIO_VALID_LICENSE) return 0;
	return m_comp.c_str();
}

//-----------------------------------------------------------------------------
void LoadLicenseFile()
{
#ifdef FEBIO_LICENSE
	FEBioLicenseKey* pkey = FEBioLicenseKey::GetInstance();
	pkey->Load();
#endif
}

//-----------------------------------------------------------------------------
// this is the function users will call
int GetLicenseKeyStatus()
{
	FEBioLicenseKey* pkey = FEBioLicenseKey::GetInstance();
	return pkey->Status();
}

//-----------------------------------------------------------------------------
// get the license user
const char* GetLicenseUser()
{
	FEBioLicenseKey* pkey = FEBioLicenseKey::GetInstance();
	return pkey->User();
}

//-----------------------------------------------------------------------------
// get the license company
const char* GetLicenseCompany()
{
	FEBioLicenseKey* pkey = FEBioLicenseKey::GetInstance();
	return pkey->Company();
}
