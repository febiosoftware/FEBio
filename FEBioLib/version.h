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
#include "febiolib_api.h"

///////////////////////////////////////////////////////////////////////////////
// FEBio version numbers
// VERSION is the main version number. This number is only incremented when
// major modifications or additions where added to the code
// SUBVERSION is only incremented when minor modifications or 
// additions where added to the code.
// SUBSUBVERSION is incremented when bugs are fixed.

#define VERSION			4
#define SUBVERSION		1
#define SUBSUBVERSION	0

///////////////////////////////////////////////////////////////////////////////
// Restart file version
// This is the version number of the restart dump file format.
// It is incremented when the structure of this file is modified.
//

#define RSTRTVERSION		0x06

namespace febio
{

FEBIOLIB_API char* getVersionString();

}
