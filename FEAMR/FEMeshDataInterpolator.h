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
#include <vector>
#include <functional>
#include <FECore/matrix.h>

// Base classes for classes that interpolate data on a mesh
class FEMeshDataInterpolator
{
public:
	FEMeshDataInterpolator();
	virtual ~FEMeshDataInterpolator();

	// do one-time initalization
	virtual bool Init();

	// set the next target point
	// should return false if the target point cannot be evaluated
	virtual bool SetTargetPoint(const vec3d& r) = 0;

	//! map source data onto target data
	//! output: tval - values at the target points
	//! input: sval - values of the source points
	virtual bool Map(std::vector<double>& tval, std::function<double(int sourceNode)> src) = 0;

	virtual double Map(int inode, std::function<double(int sourceNode)> src) = 0;
	virtual vec3d MapVec3d(int inode, std::function<vec3d(int sourceNode)> src) = 0;

	double Map(std::function<double(int sourceNode)> f);
	vec3d MapVec3d(std::function<vec3d(int sourceNode)> f);
};

inline double FEMeshDataInterpolator::Map(std::function<double(int sourceNode)> f) { return Map(0, f); }
inline vec3d FEMeshDataInterpolator::MapVec3d(std::function<vec3d(int sourceNode)> f) { return MapVec3d(0, f); }

