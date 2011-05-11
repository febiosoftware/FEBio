#pragma once

//-------------------------------------------------------------------
// Call this function to load the license file. The status of the
// license file can then be tested with the GetLicenseKeyStatus function
void LoadLicenseFile();

//-------------------------------------------------------------------
// use this function to obtain the status of the license key
// return values:
// 0 = no licences file found
// 1 = valid license file
// 2 = invalid license file
//
int GetLicenseKeyStatus();

//-------------------------------------------------------------------
// get the user name
// returns 0 if the license is invalid
const char* GetLicenseUser();

//-------------------------------------------------------------------
// get the company name
// returns 0 if the license is invalid
const char* GetLicenseCompany();
