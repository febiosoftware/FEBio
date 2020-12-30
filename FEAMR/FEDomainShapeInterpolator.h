/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FEMeshDataInterpolator.h"

class FEDomain;
class FESolidElement;
class FEMesh;
class FEOctreeSearch;

//! Maps data by using element shape functions
class FEDomainShapeInterpolator : public FEMeshDataInterpolator
{
	struct Data
	{
		FESolidElement*		el;
		double				r[3];
	};

public:
	FEDomainShapeInterpolator(FEDomain* domain);
	~FEDomainShapeInterpolator();

	bool Init();

	void SetTargetPoints(const vector<vec3d>& trgPoints);
	bool SetTargetPoint(const vec3d& r) override;

	bool Map(std::vector<double>& tval, function<double(int sourceNode)> src) override;

	double Map(int inode, function<double(int sourceNode)> src) override;
	vec3d MapVec3d(int inode, function<vec3d(int sourceNode)> src) override;

private:
	FEDomain*	m_dom;
	FEMesh*	m_mesh;

	FEOctreeSearch*	m_os;
	vector<vec3d>	m_trgPoints;
	vector<Data>	m_data;
};

