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



#pragma once

//-----------------------------------------------------------------------------
// This file defines the SDK versioning. It follows the versioning of FEBio, but
// will only be modified if the SDK is no longer compatible with older versions
// of FEBio. In that case, plugins need to be recompiled to be usable with the 
// newer version of FEBio.
#define FE_SDK_MAJOR_VERSION	4
#define FE_SDK_SUB_VERSION		1
#define FE_SDK_SUBSUB_VERSION	0

//-----------------------------------------------------------------------------
// This macro needs to be exported by all plugins in the GetSDKVersion() function.
#define FE_SDK_VERSION ((FE_SDK_MAJOR_VERSION << 16) | (FE_SDK_SUB_VERSION << 8) | (FE_SDK_SUBSUB_VERSION))
