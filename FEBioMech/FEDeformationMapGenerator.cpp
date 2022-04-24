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
#include "stdafx.h"
#include "FEDeformationMapGenerator.h"
#include <FECore/FEModel.h>
#include <FECore/FENodeDataMap.h>
#include <FECore/FEDomainMap.h>

double defgrad(FESolidElement &el, std::vector<vec3d>& X, std::vector<vec3d>& u, mat3d &F, int n);

BEGIN_FECORE_CLASS(FEDeformationMapGenerator, FEElemDataGenerator)
	ADD_PARAMETER(m_nodeDisplacementMap, "node_displacement_map");
END_FECORE_CLASS();


FEDeformationMapGenerator::FEDeformationMapGenerator(FEModel* fem) : FEElemDataGenerator(fem)
{
	m_nodeMap = nullptr;
}

FEDeformationMapGenerator::~FEDeformationMapGenerator()
{

}

bool FEDeformationMapGenerator::Init()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	FEDataMap* map = mesh.FindDataMap(m_nodeDisplacementMap);
	if (map == nullptr) return false;

	m_nodeMap = dynamic_cast<FENodeDataMap*>(map);
	if (m_nodeMap == nullptr) return false;

	if (m_nodeMap->DataType() != FE_VEC3D) return false;

	return FEElemDataGenerator::Init();
}

// generate the data array for the given element set
bool FEDeformationMapGenerator::Generate(FEDomainMap& map)
{
	const FEElementSet& set = *map.GetElementSet();

	FEMesh& mesh = *set.GetMesh();

	FEDataType dataType = map.DataType();
	if (dataType != FE_MAT3D) return false;

	int storageFormat = map.StorageFormat();
	if (storageFormat != FMT_MATPOINTS) return false;

	int N = set.Elements();
	for (int i = 0; i < N; ++i)
	{
		FESolidElement* pel = dynamic_cast<FESolidElement*>(mesh.FindElementFromID(set[i]));
		if (pel == nullptr) return false;
		FESolidElement& el = *pel;

		int ne = el.Nodes();
		vector<vec3d> u(ne), X(ne);
		for (int j = 0; j < ne; j++)
		{
			u[j] = m_nodeMap->get<vec3d>(el.m_node[j]);
			X[j] = mesh.Node(el.m_node[j]).m_r0 - u[j];
		}

		int ni = el.GaussPoints();
		for (int j = 0; j < ni; ++j)
		{
			mat3d F;
			defgrad(el, X, u, F, j);
			map.setValue(i, j, F);
		}
	}
	return true;
}

double invjac0(FESolidElement &el, std::vector<vec3d>& r0, double Ji[3][3], int n)
{
	// calculate Jacobian
	double J[3][3] = { 0 };
	int neln = el.Nodes();
	for (int i = 0; i < neln; ++i)
	{
		const double& Gri = el.Gr(n)[i];
		const double& Gsi = el.Gs(n)[i];
		const double& Gti = el.Gt(n)[i];

		const double& x = r0[i].x;
		const double& y = r0[i].y;
		const double& z = r0[i].z;

		J[0][0] += Gri * x; J[0][1] += Gsi * x; J[0][2] += Gti * x;
		J[1][0] += Gri * y; J[1][1] += Gsi * y; J[1][2] += Gti * y;
		J[2][0] += Gri * z; J[2][1] += Gsi * z; J[2][2] += Gti * z;
	}

	// calculate the determinant
	double det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
		+ J[0][1] * (J[1][2] * J[2][0] - J[2][2] * J[1][0])
		+ J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

	if (det != 0.0)
	{
		// calculate the inverse jacobian
		double deti = 1.0 / det;

		Ji[0][0] = deti * (J[1][1] * J[2][2] - J[1][2] * J[2][1]);
		Ji[1][0] = deti * (J[1][2] * J[2][0] - J[1][0] * J[2][2]);
		Ji[2][0] = deti * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

		Ji[0][1] = deti * (J[0][2] * J[2][1] - J[0][1] * J[2][2]);
		Ji[1][1] = deti * (J[0][0] * J[2][2] - J[0][2] * J[2][0]);
		Ji[2][1] = deti * (J[0][1] * J[2][0] - J[0][0] * J[2][1]);

		Ji[0][2] = deti * (J[0][1] * J[1][2] - J[1][1] * J[0][2]);
		Ji[1][2] = deti * (J[0][2] * J[1][0] - J[0][0] * J[1][2]);
		Ji[2][2] = deti * (J[0][0] * J[1][1] - J[0][1] * J[1][0]);
	}

	return det;
}

double defgrad(FESolidElement &el, std::vector<vec3d>& X, std::vector<vec3d>& u, mat3d &F, int n)
{
	// calculate inverse jacobian
    double Ji[3][3];
    invjac0(el, X, Ji, n);

	// shape function derivatives
	double *Grn = el.Gr(n);
	double *Gsn = el.Gs(n);
	double *Gtn = el.Gt(n);

	// calculate deformation gradient
	F[0][0] = F[0][1] = F[0][2] = 0;
	F[1][0] = F[1][1] = F[1][2] = 0;
	F[2][0] = F[2][1] = F[2][2] = 0;
	int neln = el.Nodes();
	for (int i = 0; i < neln; ++i)
	{
		double Gri = Grn[i];
		double Gsi = Gsn[i];
		double Gti = Gtn[i];

		double x = u[i].x;
		double y = u[i].y;
		double z = u[i].z;

		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		double GX = Ji[0][0] * Gri + Ji[1][0] * Gsi + Ji[2][0] * Gti;
		double GY = Ji[0][1] * Gri + Ji[1][1] * Gsi + Ji[2][1] * Gti;
		double GZ = Ji[0][2] * Gri + Ji[1][2] * Gsi + Ji[2][2] * Gti;

		// calculate deformation gradient F
		F[0][0] += GX * x; F[0][1] += GY * x; F[0][2] += GZ * x;
		F[1][0] += GX * y; F[1][1] += GY * y; F[1][2] += GZ * y;
		F[2][0] += GX * z; F[2][1] += GY * z; F[2][2] += GZ * z;
	}

	F[0][0] += 1.0;
	F[1][1] += 1.0;
	F[2][2] += 1.0;

	double D = F.det();

	return D;
}

