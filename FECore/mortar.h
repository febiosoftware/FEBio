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
#include "FESurface.h"
#include <vector>

//-----------------------------------------------------------------------------
// Structure for representing 2D points
struct POINT2D
{
	double x, y;
};

//-----------------------------------------------------------------------------
// Calculates the intersection between two convex polygons
FECORE_API int ConvexIntersect(POINT2D* P, int n, POINT2D* Q, int m, POINT2D* R);

//-----------------------------------------------------------------------------
class FECORE_API Patch
{
public:
	struct FACET
	{
		FACET(){}
		FACET(vec3d x[3]) { r[0] = x[0]; r[1] = x[1]; r[2] = x[2]; }
		vec3d	r[3];	//!< position of nodes

		// Evaluate the spatial position using iso-parametric coordinates
		vec3d Position(double rp, double sp)
		{
			return (r[0]*(1.0 - rp - sp) + r[1]*rp + r[2]*sp);
		}

		//! calculate the area of the patch
		double Area()
		{
			vec3d n = (r[1] - r[0])^(r[2] - r[0]);
			return n.norm()*0.5;
		}
	};

public:
	Patch(int k, int l) : m_primary_facet_id(k), m_secondary_facet_id(l) {}

	Patch(const Patch& p) { 
		m_tri = p.m_tri; 
		m_primary_facet_id = p.m_primary_facet_id;
		m_secondary_facet_id = p.m_secondary_facet_id;
	}

	Patch& operator = (const Patch& p) { 
		m_tri = p.m_tri; 
		m_primary_facet_id = p.m_primary_facet_id;
		m_secondary_facet_id = p.m_secondary_facet_id;
		return *this; 
	}

	int GetPrimaryFacetID() const { return m_primary_facet_id; }
	int GetSecondaryFacetID() const { return m_secondary_facet_id; }

public:
	//! Clear the patch
	void Clear() { m_tri.clear(); }

	//! Add a facet to the patch
	void Add(vec3d x[3]) { m_tri.push_back(FACET(x)); }

	//! retrieve a facet
	FACET& Facet(int i) { return m_tri[i]; }

	//! facet count
	int Size() { return (int) m_tri.size(); }

	//! See if the facet is empty
	bool Empty() { return m_tri.empty(); }

private:
	int		m_primary_facet_id;		//!< index of primary facet
	int		m_secondary_facet_id;	//!< index of secondary facet

	std::vector<FACET>	m_tri;	//!< triangular patches
};

//-----------------------------------------------------------------------------
class FECORE_API MortarSurface
{
public:
	MortarSurface(){}

	int Patches() { return (int) m_patch.size(); }

	Patch& GetPatch(int i) { return m_patch[i]; }

	void AddPatch(const Patch& p) { m_patch.push_back(p); }

	void Clear() { m_patch.clear(); }

private:
	std::vector<Patch>	m_patch;
};

//-----------------------------------------------------------------------------
// Calculates the intersection between two segments and adds it to the patch
FECORE_API bool CalculateMortarIntersection(FESurface& ss, FESurface& ms, int k, int l, Patch& patch);

//-----------------------------------------------------------------------------
// Calculates the mortar intersection between two surfaces
FECORE_API void CalculateMortarSurface(FESurface& ss, FESurface& ms, MortarSurface& s);

//-----------------------------------------------------------------------------
// Stores the mortar surface in STL format
FECORE_API bool ExportMortar(MortarSurface& mortar, const char* szfile);
