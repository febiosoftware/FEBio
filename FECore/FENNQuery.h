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
#include "vec3d.h"
#include <vector>
#include "fecore_api.h"

class FESurface;

//-----------------------------------------------------------------------------
//! This class is a helper class to locate the nearest neighbour on a surface

class FECORE_API FENNQuery
{
public:
	struct NODE
	{
		int		i;	// index of node
		vec3d	r;	// position of node
		double	d1;	// distance to pivot 1
		double	d2;	// distance to pivot 2
	};

public:
	FENNQuery(FESurface* ps = 0);
	virtual ~FENNQuery();

	//! initialize search structures
	void Init();
	void InitReference();

	//! attach to a surface
	void Attach(FESurface* ps) { m_ps = ps; }

	//! find the neirest neighbour of r
	int Find(vec3d x);	
	int FindReference(vec3d x);	

protected:
	int FindRadius(double r);

protected:
	FESurface*	m_ps;	//!< the surface to search
	std::vector<NODE>	m_bk;	// BK tree

	vec3d	m_q1;	// pivot 1
	vec3d	m_q2;	// pivot 2

	int		m_imin;	// last found index
};

// function for finding the k closest neighbors
int FECORE_API findNeirestNeighbors(const std::vector<vec3d>& point, const vec3d& x, int k, std::vector<int>& closestNodes);
