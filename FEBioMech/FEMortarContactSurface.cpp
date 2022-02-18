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
#include "FEMortarContactSurface.h"

//-----------------------------------------------------------------------------
FEMortarContactSurface::FEMortarContactSurface(FEModel* pfem) : FEContactSurface(pfem)
{
}

//-----------------------------------------------------------------------------
bool FEMortarContactSurface::Init()
{
	if (FEContactSurface::Init() == false) return false;

	int NN = Nodes();
	m_gap.resize(NN, vec3d(0,0,0));
	return true;
}

//-----------------------------------------------------------------------------
void FEMortarContactSurface::UpdateNodalAreas()
{
	int NN = Nodes();
	int NF = Elements();
	m_A.resize(NN, 0.0);

	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = Element(i);
		double a = FaceArea(el);

		int nn = el.Nodes();
		double fa = a / (double) nn;
		for (int j=0; j<nn; ++j) m_A[el.m_lnode[j]] += fa;
	}

	for (int i=0; i<NN; ++i) m_A[i] = 1.0/m_A[i];
}
