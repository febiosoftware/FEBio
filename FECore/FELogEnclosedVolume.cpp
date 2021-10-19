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
#include "FELogEnclosedVolume.h"
#include "FESurface.h"
#include "FENode.h"

BEGIN_FECORE_CLASS(FELogEnclosedVolume, FELogSurfaceData)
END_FECORE_CLASS();

double FELogEnclosedVolume::value(FESurface& surface)
{
	// loop over all elements
	double vol = 0.0;
	int NE = surface.Elements();
	vec3d x[FEElement::MAX_NODES];
	for (int i = 0; i < NE; ++i)
	{
		// get the next element
		FESurfaceElement& el = surface.Element(i);

		// get the nodal coordinates
		int neln = el.Nodes();
		for (int j = 0; j < neln; ++j) x[j] = surface.Node(el.m_lnode[j]).m_rt;

		// loop over integration points
		double* w = el.GaussWeights();
		int nint = el.GaussPoints();
		for (int n = 0; n < nint; ++n)
		{
			// evaluate the position vector at this point
			vec3d r = el.eval(x, n);

			// calculate the tangent vectors
			double* Gr = el.Gr(n);
			double* Gs = el.Gs(n);
			vec3d dxr(0, 0, 0), dxs(0, 0, 0);
			for (int j = 0; j < neln; ++j)
			{
				dxr += x[j] * Gr[j];
				dxs += x[j] * Gs[j];
			}

			// update volume
			vol += w[n] * (r * (dxr ^ dxs));
		}
	}
	return vol / 3.0;
}
