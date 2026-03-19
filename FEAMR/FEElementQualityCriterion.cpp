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
#include "FEElementQualityCriterion.h"
#include <FECore/FEElement.h>
#include <FECore/FEMeshPartition.h>
#include <FECore/FEMesh.h>

BEGIN_FECORE_CLASS(FEMeanRatioQualityCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(minQuality, "min_quality")->setLongName("Minimum element quality");
END_FECORE_CLASS();

FEMeanRatioQualityCriterion::FEMeanRatioQualityCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	minQuality = 1.0;
}

bool FEMeanRatioQualityCriterion::GetElementValue(FEElement& el, double& value)
{
	if ((el.Shape() != ET_TET4) && (el.Shape() != ET_HEX8) && (el.Shape() != ET_PENTA6)) return false;
	value = ElementQuality(el);
	return (value < minQuality);
}

// Edge list for different element types
static const int hexCornerEdges[8][3][2] =
{
	{ {0,1}, { 0, 3}, {0,4} },
	{ {1,2}, { 1,0 }, {1,5} },
	{ {2,3}, { 2,1 }, {2,6} },
	{ {3,0}, { 3,2 }, {3,7} },
	{ {4,0}, { 4,7 }, {4,5} },
	{ {5,1}, { 5,4 }, {5,6} },
	{ {6,2}, { 6,5 }, {6,7} },
	{ {7,3}, { 7,6 }, {7,4} } 
};

static const int wedgeCornerEdges[6][3][2] =
{
	{ {0,1}, {0,2}, {0,3} },
	{ {1,2}, {1,0}, {1,4} },
	{ {2,0}, {2,1}, {2,5} },
	{ {3,4}, {3,0}, {3,5} },
	{ {4,5}, {4,1}, {4,3} },
	{ {5,3}, {5,2}, {5,4} }
};

static const int tetCornerEdges[4][3][2] =
{
	{ {0,1}, {0,2}, {0,3} },
	{ {1,2}, {1,0}, {1,3} },
	{ {2,0}, {2,1}, {2,3} },
	{ {3,0}, {3,2}, {3,1} }
};

double FEMeanRatioQualityCriterion::ElementQuality(FEElement& el)
{
	FEMeshPartition* partition = el.GetMeshPartition();
	if (partition == nullptr) return 0.0;

	FEMesh& mesh = *partition->GetMesh();

	// get the nodal coordinates
	vec3d r[FEElement::MAX_NODES];
	for (int i = 0; i < el.Nodes(); ++i) r[i] = mesh.Node(el.m_node[i]).m_rt;

	// for each node, find 3 edges that meet at the node
	int n = el.Nodes();
	const int(*edgeData)[3][2] = nullptr;
	switch (el.Shape())
	{
	case ET_HEX8  : edgeData = hexCornerEdges; break;
	case ET_PENTA6: edgeData = wedgeCornerEdges; break;
	case ET_TET4  : edgeData = tetCornerEdges; break;
	default: return 1.0;
	}

	// calculate mean ratio at each node and store minimum value
	double minMR = 1e99;
	vec3d e[3];
	double L[3];
	for (int i = 0; i < n; ++i)
	{
		// get the 3 edges that meet at this node
		for (int k = 0; k < 3; ++k)
		{
			int a = edgeData[i][k][0];
			int b = edgeData[i][k][1];
			vec3d va = r[a];
			vec3d vb = r[b];
			e[k] = vb - va;
			L[k] = e[k].Length();
		}

		// get the jacobian determinant at this node
		double detJ = e[0] * (e[1] ^ e[2]);
		double JF2 = e[0].SqrLength() + e[1].SqrLength() + e[2].SqrLength();

		// calculate the mean ratio
		double MR = 3.0 * pow(detJ, 2.0 / 3.0) / JF2;
		if (MR < minMR) minMR = MR;
	}

	return minMR;
}

//==================================================================

BEGIN_FECORE_CLASS(FEScaledJacobianQualityCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(minQuality, "min_quality")->setLongName("Minimum element quality");
END_FECORE_CLASS();

FEScaledJacobianQualityCriterion::FEScaledJacobianQualityCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	minQuality = 1.0;
}

bool FEScaledJacobianQualityCriterion::GetElementValue(FEElement& el, double& value)
{
	if ((el.Shape() != ET_TET4) && (el.Shape() != ET_HEX8) && (el.Shape() != ET_PENTA6)) return false;
	value = ElementQuality(el);
	return (value < minQuality);
}

double FEScaledJacobianQualityCriterion::ElementQuality(FEElement& el)
{
	FEMeshPartition* partition = el.GetMeshPartition();
	if (partition == nullptr) return 0.0;

	FEMesh& mesh = *partition->GetMesh();

	// get the nodal coordinates
	vec3d r[FEElement::MAX_NODES];
	for (int i = 0; i < el.Nodes(); ++i) r[i] = mesh.Node(el.m_node[i]).m_rt;

	// for each node, find 3 edges that meet at the node
	int n = el.Nodes();
	const int(*edgeData)[3][2] = nullptr;
	switch (el.Type())
	{
	case ET_HEX8  : edgeData = hexCornerEdges; break;
	case ET_PENTA6: edgeData = wedgeCornerEdges; break;
	case ET_TET4  : edgeData = tetCornerEdges; break;
	default: return 0.0;
	}

	// calculate scaled jacobian at each node and store minimum value
	double minSJ = 1e99;
	vec3d e[3];
	double L[3];
	for (int i = 0; i < n; ++i)
	{
		// get the 3 edges that meet at this node
		for (int k = 0; k < 3; ++k)
		{
			int a = edgeData[i][k][0];
			int b = edgeData[i][k][1];
			vec3d va = r[a];
			vec3d vb = r[b];
			e[k] = vb - va;
			L[k] = e[k].Length();
		}

		// get the jacobian determinant at this node
		double detJ = e[0] * (e[1] ^ e[2]);

		// calculate the scaled jacobian
		double SJ = detJ / (L[0] * L[1] * L[2]);
		if (SJ < minSJ) minSJ = SJ;
	}

	return minSJ;
}