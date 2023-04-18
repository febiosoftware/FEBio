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
#include "FEDataExport.h"
#include "vec3d.h"
using namespace std;

//-----------------------------------------------------------------------------
void FEDataExport::Serialize(FEDataStream& ar)
{
	if ((m_type == PLT_VEC3F)&&(m_fmt == FMT_NODE))
	{
		vector<vec3d>& v = *(static_cast<vector<vec3d>*>(m_pd));

		int n = (int) v.size();
		for (int i=0; i<n; ++i) ar << v[i];
	}
	else if ((m_type == PLT_FLOAT)&&(m_fmt == FMT_REGION))
	{
		double& d = *(static_cast<double*>(m_pd));
		ar << d;
	}
	else if ((m_type == PLT_FLOAT) && (m_fmt == FMT_NODE))
	{
		vector<double>& v = *(static_cast<vector<double>*>(m_pd));

		int n = (int)v.size();
		for (int i = 0; i<n; ++i) ar << v[i];
	}
	else
	{
		assert(false);
	}
}
