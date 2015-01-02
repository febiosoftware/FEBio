#pragma once

//-----------------------------------------------------------------------------
// This file defines the SDK versioning. It follows the versioning of FEBio, but
// will only be modified if the SDK is no longer compatible with older versions
// of FEBio. In that case, plugins need to be recompiled to be usable with the 
// newer version of FEBio.
#define FE_SDK_MAJOR_VERSION	2
#define FE_SDK_SUB_VERSION		3
#define FE_SDK_SUBSUB_VERSION	0

//-----------------------------------------------------------------------------
// This macro needs to be exported by all plugins in the GetSDKVersion() function.
#define FE_SDK_VERSION ((FE_SDK_MAJOR_VERSION << 16) | (FE_SDK_SUB_VERSION << 8) | (FE_SDK_SUBSUB_VERSION))
