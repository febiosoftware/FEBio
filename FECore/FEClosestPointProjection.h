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
#include "FENNQuery.h"
#include "FEElemElemList.h"
#include "FENodeElemList.h"

//-----------------------------------------------------------------------------
// This class can be used to find the closest point projection of a point
// onto a surface.
class FECORE_API FEClosestPointProjection
{
public:
	//! constructor
	FEClosestPointProjection(FESurface& s);

	//! Initialization
	bool Init();

	//! Project a point onto surface
	FESurfaceElement* Project(const vec3d& x, vec3d& q, vec2d& r);

	//! Project a node onto a surface
	FESurfaceElement* Project(int n, vec3d& q, vec2d& r);

	//! Project a point of a surface element onto a surface
	FESurfaceElement* Project(FESurfaceElement* pse, int intgrPoint, vec3d& q, vec2d& r);

public:
	//! Set the projection tolerance
	void SetTolerance(double t) { m_tol = t; }

	//! get the projection tolerance
	double GetTolerance() { return m_tol; }

	//! Set the search radius (used in self-projection)
	void SetSearchRadius(double s) { m_rad = s; }

	//! set if the projection should handle special cases
	void HandleSpecialCases(bool b) { m_bspecial = b; }

	//! set if boundary projections are allowed
	void AllowBoundaryProjections(bool b) { m_projectBoundary = b; }

private:
	bool ContainsElement(FESurfaceElement* el);
	FESurfaceElement* ProjectSpecial(int closestPoint, const vec3d& x, vec3d& q, vec2d& r);

protected:
	double	m_tol;	//!< projection tolerance
	double	m_rad;	//!< search radius
	bool	m_bspecial;	//!< try to handle special cases
	bool	m_projectBoundary;	//!< allow boundary projections

protected:
	FESurface&		m_surf;		//!< reference to surface
	FENNQuery		m_SNQ;		//!< used to find the nearest neighbour
	FENodeElemList	m_NEL;		//!< node-element tree
	FEElemElemList	m_EEL;		//!< element neighbor list
};
