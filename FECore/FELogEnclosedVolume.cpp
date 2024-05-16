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
	double vol = 0.0;
	int NE = surface.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FESurfaceElement& el = surface.Element(i);
		double* w = el.GaussWeights();
		for (int n = 0; n < el.GaussPoints(); ++n)
		{
			vec3d xi = surface.Local2Global(el, n);
			vec3d g[2];
			surface.CoBaseVectors(el, n, g);
			double DVj = xi * (g[0] ^ g[1]) / 3;
			vol += DVj * w[n];
		}
	}
	return vol;
}

double FELogEnclosedVolumeChange::value(FESurface& surface)
{
	double DV = 0.0;
	int NE = surface.Elements();
	for (int i = 0; i < NE; ++i) {
		FESurfaceElement& el = surface.Element(i);
		double* w = el.GaussWeights();
		for (int n = 0; n < el.GaussPoints(); ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);
			vec3d xi = surface.Local2Global(el, n);
			vec3d g[2];
			surface.CoBaseVectors(el, n, g);
			vec3d Xi = surface.Local2Global0(el, n);
			vec3d G[2];
			surface.CoBaseVectors0(el, n, G);
			double DVj = (xi * (g[0] ^ g[1]) - Xi * (G[0] ^ G[1])) / 3;
			DV += DVj * w[n];
		}
	}
	return DV;
}
