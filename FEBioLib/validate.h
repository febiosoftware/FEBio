#pragma once
#include <string>
#include "febiolib_api.h"

//-------------------------------------------------------------------
// use this function to obtain the status of the license key
// return values:
// 0 = non-commercial license
// 1 = valid commercial license file
// 2 = invalid license file
//
FEBIOLIB_API int GetLicenseKeyStatus(const char* licenseKey = 0);

//-------------------------------------------------------------------
// load the license key 
// This will try to find a file named license.txt in the application path
// and read whatever is inside that file. 
// returns an empty string if the license file was not found
FEBIOLIB_API std::string LoadLicenseKey();
