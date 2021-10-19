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
#include "FESurface.h"
#include "FEOctree.h"

//-----------------------------------------------------------------------------
//! This class calculates the normal projection on to a surface.
//! This is used by some contact algorithms.
class FECORE_API FENormalProjection
{
public:
	//! constructor
	FENormalProjection(FESurface& s);

	// initialization
	void Init();

	void SetTolerance(double tol) { m_tol = tol; }
	void SetSearchRadius(double srad) { m_rad = srad; }

public:
	//! find the intersection of a ray with the surface
	FESurfaceElement* Project(vec3d r, vec3d n, double rs[2]);
	FESurfaceElement* Project2(vec3d r, vec3d n, double rs[2]);
	FESurfaceElement* Project3(const vec3d& r, const vec3d& n, double rs[2], int* pei = 0);

	vec3d Project(const vec3d& r, const vec3d& N);
	vec3d Project2(const vec3d& r, const vec3d& N);

private:
	double	m_tol;	//!< projection tolerance
	double	m_rad;	//!< search radius

private:
	FESurface&	m_surf;	//!< the target surface
	FEOctree	m_OT;	//!< used to optimize ray-surface intersections
};
