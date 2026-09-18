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
#include "FELogData.h"
#include <assert.h>

FELogData::FELogData(FEModel* fem) : FECoreBase(fem)
{

}

int Vec3dLogTraits::ComponentIndex(const std::string& c)
{
	if (c == "x") return 0;
	if (c == "y") return 1;
	if (c == "z") return 2;
	assert(false);
	return -1;
}

double Vec3dLogTraits::Component(const vec3d& v, int i)
{
	switch (i)
	{
	case 0: return v.x;
	case 1: return v.y;
	case 2: return v.z;
	default:
		assert(false);
		return 0.0;
	}
}

int Mat3dsLogTraits::ComponentIndex(const std::string& c)
{
	if (c == "xx") return 0;
	if (c == "yy") return 1;
	if (c == "zz") return 2;
	if (c == "xy" || c == "yx") return 3;
	if (c == "yz" || c == "zy") return 4;
	if (c == "xz" || c == "zx") return 5;
	assert(false);
	return -1;
}

double Mat3dsLogTraits::Component(const mat3ds& v, int i)
{
	switch (i)
	{
	case 0: return v.xx();
	case 1: return v.yy();
	case 2: return v.zz();
	case 3: return v.xy();
	case 4: return v.yz();
	case 5: return v.xz();
	default:
		assert(false);
		return 0.0;
	}
}

int Mat3dLogTraits::ComponentIndex(const std::string& c)
{
	if (c == "xx") return 0;
	if (c == "xy") return 1;
	if (c == "xz") return 2;
	if (c == "yx") return 3;
	if (c == "yy") return 4;
	if (c == "yz") return 5;
	if (c == "zx") return 6;
	if (c == "zy") return 7;
	if (c == "zz") return 8;
	assert(false);
	return -1;
}

double Mat3dLogTraits::Component(const mat3d& v, int i)
{
	switch (i)
	{
	case 0: return v(0,0);
	case 1: return v(0,1);
	case 2: return v(0,2);
	case 3: return v(1,0);
	case 4: return v(1,1);
	case 5: return v(1,2);
	case 6: return v(2,0);
	case 7: return v(2,1);
	case 8: return v(2,2);
	default:
		assert(false);
		return 0.0;
	}
}
