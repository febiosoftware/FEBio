/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "writeplot.h"
#include "FESPRProjection.h"

//-------------------------------------------------------------------------------------------------
void writeSPRElementValueMat3dd(FESolidDomain& dom, FEDataStream& ar, std::function<mat3dd(const FEMaterialPoint&)> fnc, int interpolOrder)
{
	int NN = dom.Nodes();
	int NE = dom.Elements();

	// build the element data array
	vector< vector<double> > ED[3];
	ED[0].resize(NE);
	ED[1].resize(NE);
	ED[2].resize(NE);
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& e = dom.Element(i);
		int nint = e.GaussPoints();
		ED[0][i].assign(nint, 0.0);
		ED[1][i].assign(nint, 0.0);
		ED[2][i].assign(nint, 0.0);
	}

	// this array will store the results
	FESPRProjection map;
	map.SetInterpolationOrder(interpolOrder);
	vector<double> val[3];

	// fill the ED array
	for (int i = 0; i < NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		for (int j = 0; j < nint; ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			mat3dd v = fnc(mp);

			ED[0][i][j] = v.diag(0);
			ED[1][i][j] = v.diag(1);
			ED[2][i][j] = v.diag(2);
		}
	}

	// project to nodes
	map.Project(dom, ED[0], val[0]);
	map.Project(dom, ED[1], val[1]);
	map.Project(dom, ED[2], val[2]);

	// copy results to archive
	for (int i = 0; i<NN; ++i)
	{
		ar.push_back((float)val[0][i]);
		ar.push_back((float)val[1][i]);
		ar.push_back((float)val[2][i]);
	}
}

//-------------------------------------------------------------------------------------------------
void writeSPRElementValueMat3ds(FESolidDomain& dom, FEDataStream& ar, std::function<mat3ds(const FEMaterialPoint&)> fnc, int interpolOrder)
{
	const int LUT[6][2] = { { 0,0 },{ 1,1 },{ 2,2 },{ 0,1 },{ 1,2 },{ 0,2 } };

	int NN = dom.Nodes();
	int NE = dom.Elements();

	// build the element data array
	vector< vector<double> > ED[6];
	for (int n = 0; n < 6; ++n)
	{
		ED[n].resize(NE);
		for (int i = 0; i < NE; ++i)
		{
			FESolidElement& e = dom.Element(i);
			int nint = e.GaussPoints();
			ED[n][i].assign(nint, 0.0);
		}
	}

	// this array will store the results
	FESPRProjection map;
	map.SetInterpolationOrder(interpolOrder);
	vector<double> val[6];

	// fill the ED array
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			mat3ds s = fnc(mp);

			// loop over stress components
			for (int n = 0; n < 6; ++n)
			{
				ED[n][i][j] = s(LUT[n][0], LUT[n][1]);
			}
		}
	}

	// project to nodes
	// loop over stress components
	for (int n = 0; n<6; ++n)
	{
		map.Project(dom, ED[n], val[n]);
	}

	// copy results to archive
	for (int i = 0; i<NN; ++i)
	{
		ar.push_back((float)val[0][i]);
		ar.push_back((float)val[1][i]);
		ar.push_back((float)val[2][i]);
		ar.push_back((float)val[3][i]);
		ar.push_back((float)val[4][i]);
		ar.push_back((float)val[5][i]);
	}
}
