/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEBCPrescribedDeformation.h"
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FEModel.h>
#include "FEBioMech.h"

BEGIN_FECORE_CLASS(FEBCPrescribedDeformation, FEPrescribedNodeSet)
	ADD_PARAMETER(m_scale, "scale");
	ADD_PARAMETER(m_F    , "F");
END_FECORE_CLASS();

FEBCPrescribedDeformation::FEBCPrescribedDeformation(FEModel* pfem) : FEPrescribedNodeSet(pfem)
{
	m_scale = 1.0;
	m_F.unit();
}

//-----------------------------------------------------------------------------
bool FEBCPrescribedDeformation::SetDofList(FEDofList& dofs)
{
	return dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
}

//-----------------------------------------------------------------------------
void FEBCPrescribedDeformation::CopyFrom(FEBoundaryCondition* pbc)
{
	FEBCPrescribedDeformation* ps = dynamic_cast<FEBCPrescribedDeformation*>(pbc); assert(ps);
	m_scale = ps->m_scale;
	m_F = ps->m_F;
	CopyParameterListState(ps->GetParameterList());
}

//-----------------------------------------------------------------------------
void FEBCPrescribedDeformation::SetDeformationGradient(const mat3d& F)
{
	m_F = F;
}

//-----------------------------------------------------------------------------
void FEBCPrescribedDeformation::GetNodalValues(int nodelid, std::vector<double>& val)
{
	vec3d X = GetNodeSet()->Node(nodelid)->m_r0;
	mat3ds XX = dyad(X);
	vec3d x = m_F*X;
	vec3d u = (x - X)*m_scale;
	
	val[0] = u.x;
	val[1] = u.y;
	val[2] = u.z;
}

//=============================================================================
BEGIN_FECORE_CLASS(FEBCPrescribedDeformation2O, FEPrescribedNodeSet)
	ADD_PARAMETER(m_scale, "scale");
	ADD_PARAMETER(m_F    , "F");
	ADD_PARAMETER(m_G    , "G");
	ADD_PARAMETER(m_refNode, "reference");
END_FECORE_CLASS();

FEBCPrescribedDeformation2O::FEBCPrescribedDeformation2O(FEModel* pfem) : FEPrescribedNodeSet(pfem)
{
	m_scale = 1.0;
	m_F.unit();
	m_G.zero();
	m_refNode = -1;
}

//-----------------------------------------------------------------------------
bool FEBCPrescribedDeformation2O::SetDofList(FEDofList& dofs)
{
	return dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
}

//-----------------------------------------------------------------------------
bool FEBCPrescribedDeformation2O::Init()
{
	if (m_refNode < 0) return false;
	return FEPrescribedNodeSet::Init();
}

//-----------------------------------------------------------------------------
// Sets the displacement scale factor. An optional load curve index can be given
// of the load curve that will control the scale factor.
void FEBCPrescribedDeformation2O::SetScale(double s, int lc)
{
	m_scale = s;
	if (lc >= 0)
	{
		FEParam& p = *FEParamContainer::FindParameterFromData((void*)(&m_scale));
		GetFEModel()->AttachLoadController(&p, lc);
	}
}

//-----------------------------------------------------------------------------
void FEBCPrescribedDeformation2O::CopyFrom(FEBoundaryCondition* pbc)
{
	FEBCPrescribedDeformation2O* ps = dynamic_cast<FEBCPrescribedDeformation2O*>(pbc); assert(ps);
	m_scale = ps->m_scale;
	m_F = ps->m_F;
	m_G = ps->m_G;
	m_refNode = ps->m_refNode;
	CopyParameterListState(ps->GetParameterList());
}

//-----------------------------------------------------------------------------
void FEBCPrescribedDeformation2O::SetReferenceNode(int n)
{
	m_refNode = n;
}

//-----------------------------------------------------------------------------
void FEBCPrescribedDeformation2O::SetDeformationGradient(const mat3d& F)
{
	m_F = F;
}

//-----------------------------------------------------------------------------
void FEBCPrescribedDeformation2O::SetDeformationHessian(const tens3drs& G)
{
	m_G = G;
}

//-----------------------------------------------------------------------------
void FEBCPrescribedDeformation2O::GetNodalValues(int nodelid, std::vector<double>& val)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	vec3d X1 = mesh.Node(m_refNode).m_r0;

	vec3d X = GetNodeSet()->Node(nodelid)->m_r0;

	mat3ds XX = dyad(X);
	mat3ds XX1 = dyad(X1);
	mat3d U = m_F - mat3dd(1.0);
	vec3d u = U*(X - X1) + m_G.contract2s(XX - XX1)*0.5;
	u*=m_scale;

	val[0] = u.x;
	val[1] = u.y;
	val[2] = u.z;
}
