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

//-----------------------------------------------------------------------------
//  This class stores the coordinates of a bounding box
//
class FECORE_API FEBoundingBox
{
public:
	FEBoundingBox() {}
	FEBoundingBox(const vec3d& x) : r0(x), r1(x) {}
	FEBoundingBox(const vec3d& x0, const vec3d& x1) : r0(x0), r1(x1) {}

	// center of box
	vec3d center() const { return (r0 + r1)*0.5; }

	// dimensions of box
	double width() const { return (r1.x - r0.x); }
	double height() const { return (r1.y - r0.y); }
	double depth() const { return (r1.z - r0.z); }

	// max dimension
	double radius() const
	{
		double w = width();
		double h = height();
		double d = depth();

		if ((w >= d) && (w >= h)) return w;
		if ((h >= w) && (h >= d)) return h;

		return d;
	}

	// add a point and grow the box if necessary
	void add(const vec3d& r)
	{
		if (r.x < r0.x) r0.x = r.x;
		if (r.y < r0.y) r0.y = r.y;
		if (r.z < r0.z) r0.z = r.z;
		if (r.x > r1.x) r1.x = r.x;
		if (r.y > r1.y) r1.y = r.y;
		if (r.z > r1.z) r1.z = r.z;
	}

	// inflate the box
	void inflate(double dx, double dy, double dz)
	{
		r0.x -= dx; r1.x += dx;
		r0.y -= dy; r1.y += dy;
		r0.z -= dz; r1.z += dz;
	}

	// translate the box
	void translate(const vec3d& t)
	{
		r0 += t;
		r1 += t;
	}

	// check whether a point is inside or not
	bool IsInside(const vec3d& r) const
	{
		return ((r.x >= r0.x) && (r.y >= r0.y) && (r.z >= r0.z) && (r.x <= r1.x) && (r.y <= r1.y) && (r.z <= r1.z));
	}

private:
	vec3d	r0, r1; // coordinates of opposite corners
};
