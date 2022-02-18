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

void ProjectToNodes(FEDomain& dom, vector<double>& nodeVals, function<double(FEMaterialPoint& mp)> f)
{
	// temp storage 
	double si[FEElement::MAX_INTPOINTS];
	double sn[FEElement::MAX_NODES];

	// allocate nodeVals and create valence array (tag)
	int NN = dom.Nodes();
	vector<int> tag(NN, 0);
	nodeVals.assign(NN, 0.0);

	// loop over all elements
	int NE = dom.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int ne = e.Nodes();
		int ni = e.GaussPoints();

		// get the integration point values
		for (int k = 0; k < ni; ++k)
		{
			FEMaterialPoint& mp = *e.GetMaterialPoint(k);
			double v = f(mp);
			si[k] = v;
		}

		// project to nodes
		e.project_to_nodes(si, sn);

		for (int j = 0; j < ne; ++j)
		{
			nodeVals[e.m_lnode[j]] += sn[j];
			tag[e.m_lnode[j]]++;
		}
	}

	for (int i = 0; i < NN; ++i)
	{
		if (tag[i] > 0) nodeVals[i] /= (double)tag[i];
	}
}

void writeRelativeError(FEDomain& dom, FEDataStream& a, function<double(FEMaterialPoint& mp)> f)
{
	int NE = dom.Elements();
	int NN = dom.Nodes();

	// calculate the recovered nodal values
	vector<double> sn(NN);
	ProjectToNodes(dom, sn, f);

	// find the min and max values
	double smin = 1e99, smax = -1e99;
	for (int i = 0; i < NE; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int ni = el.GaussPoints();
		for (int j = 0; j < ni; ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			double sj = f(mp);

			if (sj < smin) smin = sj;
			if (sj > smax) smax = sj;
		}
	}
	if (fabs(smin - smax) < 1e-12) smax++;

	// calculate errors
	double ev[FEElement::MAX_NODES];
	for (int i = 0; i < NE; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int ne = el.Nodes();
		int ni = el.GaussPoints();

		// get the nodal values
		for (int j = 0; j < ne; ++j)
		{
			ev[j] = sn[el.m_lnode[j]];
		}

		// evaluate element error
		double max_err = 0;
		for (int j = 0; j < ni; ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			double sj = f(mp);

			double snj = el.Evaluate(ev, j);

			double err = fabs(sj - snj) / (smax - smin);
			if (err > max_err) max_err = err;
		}

		a << max_err;
	}
}
