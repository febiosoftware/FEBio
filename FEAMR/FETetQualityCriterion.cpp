/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2026 University of Utah, The Trustees of Columbia University in
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
#include "FETetQualityCriterion.h"
#include <FECore/FEElement.h>
#include <FECore/FEMeshPartition.h>
#include <FECore/FEMesh.h>


BEGIN_FECORE_CLASS(FETetQualityCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(minQuality, "min_quality")->setLongName("Minimum tet quality");
END_FECORE_CLASS();

FETetQualityCriterion::FETetQualityCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	minQuality = 0.0;
}

bool FETetQualityCriterion::GetElementValue(FEElement& el, double& value)
{
	if (el.Shape() != ET_TET4) return false;
	value = TetQuality(el);
	return (value >= minQuality);
}

double FETetQualityCriterion::TetQuality(FEElement& el)
{
	FEMeshPartition* partition = el.GetMeshPartition();
	if (partition == nullptr) return 0.0;

	FEMesh& mesh = *partition->GetMesh();

	// get the tet's nodal coordinates
	vec3d p[4];
	for (int i = 0; i < 4; ++i) p[i] = mesh.Node(el.m_node[i]).m_rt;

	// setup system of equation
	mat3d A;
	A[0][0] = p[1].x - p[0].x; A[0][1] = p[1].y - p[0].y; A[0][2] = p[1].z - p[0].z;
	A[1][0] = p[2].x - p[0].x; A[1][1] = p[2].y - p[0].y; A[1][2] = p[2].z - p[0].z;
	A[2][0] = p[3].x - p[0].x; A[2][1] = p[3].y - p[0].y; A[2][2] = p[3].z - p[0].z;
	A = A.inverse();

	// setup RHS
	vec3d b;
	b.x = 0.5 * (p[1].x * p[1].x - p[0].x * p[0].x + p[1].y * p[1].y - p[0].y * p[0].y + p[1].z * p[1].z - p[0].z * p[0].z);
	b.y = 0.5 * (p[2].x * p[2].x - p[0].x * p[0].x + p[2].y * p[2].y - p[0].y * p[0].y + p[2].z * p[2].z - p[0].z * p[0].z);
	b.z = 0.5 * (p[3].x * p[3].x - p[0].x * p[0].x + p[3].y * p[3].y - p[0].y * p[0].y + p[3].z * p[3].z - p[0].z * p[0].z);

	// find the center of the circum sphere
	vec3d c = A * b;

	// find the radius of the circum sphere
	double R2 = (p[0].x - c.x) * (p[0].x - c.x) + (p[0].y - c.y) * (p[0].y - c.y) + (p[0].z - c.z) * (p[0].z - c.z);
	double R = sqrt(R2);

	// find the shortest edge
	const int ET[6][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 0, 3 }, { 1, 3 }, { 2, 3 } };

	double L2, L2min = 1e99;
	for (int i = 0; i < 6; ++i)
	{
		int j = ET[i][0];
		int k = ET[i][1];
		L2 = (p[j].x - p[k].x) * (p[j].x - p[k].x) + (p[j].y - p[k].y) * (p[j].y - p[k].y) + (p[j].z - p[k].z) * (p[j].z - p[k].z);
		if (L2 < L2min) L2min = L2;
	}
	double L = sqrt(L2min);

	return R / L;
}
