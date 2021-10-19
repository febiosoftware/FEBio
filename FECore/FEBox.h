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
#include "fecore_api.h"

//-----------------------------------------------------------------------------
class FEMesh;
class FEMeshPartition;

//-----------------------------------------------------------------------------
// Helper class for finding bounding boxes.
class FECORE_API FEBox
{
public:
	FEBox();
	FEBox(const vec3d& r0, const vec3d& r1);

	FEBox(const FEMesh& mesh);
	FEBox(const FEMeshPartition& dom);

	vec3d center() { return (m_r0 + m_r1)*0.5; }
	
	double width () { return m_r1.x - m_r0.x; };	//!< x-size
	double height() { return m_r1.y - m_r0.y; };	//!< y-size
	double depth () { return m_r1.z - m_r0.z; };	//!< z-size

	// the maximum size
	double maxsize();

	// union of box and node
	void add(const vec3d& r);

private:
	vec3d	m_r0, m_r1;
};
