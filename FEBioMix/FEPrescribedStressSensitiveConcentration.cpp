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
#include "FEPrescribedStressSensitiveConcentration.h"
#include <FECore/FEMesh.h>
#include <FEBioMech/FEBioMech.h>
#include "FEBioMix.h"
#include <iostream>
#include <FEBioMech/FEElasticMaterialPoint.h>
#include <FECore/FENodeElemList.h>
#include <FEBioMix/FESolutesMaterialPoint.h>
#include <FEBioMix/FEBiphasic.h>
#include <FEBioMech/FEKinematicGrowth.h>

BEGIN_FECORE_CLASS(FEPrescribedStressSensitiveConcentration, FEPrescribedDOF)
	ADD_PARAMETER(m_dof, "dof", 0, "$(dof_list:concentration)");
	ADD_PARAMETER(m_value, "value");
	ADD_PARAMETER(stress0, "residual_stress")->setLongName("initial residual stress");
	ADD_PARAMETER(m_a0, "min_stress_scale");
	ADD_PARAMETER(m_a, "stress_scale_amplitude");
	ADD_PARAMETER(m_b, "stress_scale_width");
END_FECORE_CLASS();

FEPrescribedStressSensitiveConcentration::FEPrescribedStressSensitiveConcentration(FEModel* pfem) : FEPrescribedDOF(pfem)
{

}

bool FEPrescribedStressSensitiveConcentration::Init()
{
	if (FEPrescribedDOF::Init() == false) return false;
	return true;
}

//-----------------------------------------------------------------------------
mat3ds FEPrescribedStressSensitiveConcentration::GetStress(FEElement& m_elem, int node_id)
{
	mat3ds si;
	for (int i = 0; i < m_elem.GaussPoints(); ++i) {
		FEMaterialPoint* pt = m_elem.GetMaterialPoint(i);
		FEElasticMaterialPoint* ep = pt->ExtractData<FEElasticMaterialPoint>();
		si += ep->m_s;
	}
	si /= m_elem.GaussPoints();
	return si;
}

//-----------------------------------------------------------------------------
double FEPrescribedStressSensitiveConcentration::GetConcentration(FEElement& m_elem, int node_id)
{
	double ci = 0.0;
	for (int i = 0; i < m_elem.GaussPoints(); ++i) {
		FEMaterialPoint& pt = *m_elem.GetMaterialPoint(i);
		ci += m_value(pt);
	}
	ci /= m_elem.GaussPoints();
	return ci;
}

//-----------------------------------------------------------------------------
mat3ds FEPrescribedStressSensitiveConcentration::GetNodalStress(int node_id)
{
	// get the mesh to which this surface belongs
	FENodeElemList& NEL = GetMesh().NodeElementList();
	int nval = NEL.Valence(node_id);
	FEElement** ppe = NEL.ElementList(node_id);

	// Get the elements that the node belongs to
	mat3ds s = mat3ds(0.0);
	for (int i = 0; i < nval; ++i)
	{
		s += GetStress(*ppe[i], node_id);
	}
	s = s / nval;
	return s;
}

//-----------------------------------------------------------------------------
double FEPrescribedStressSensitiveConcentration::GetNodalConcentration(int node_id)
{
	// get the mesh to which this surface belongs
	FENodeElemList& NEL = GetMesh().NodeElementList();
	int nval = NEL.Valence(node_id);
	FEElement** ppe = NEL.ElementList(node_id);

	// Get the elements that the node belongs to
	double s = 0.0;
	//double vol = 0.0;
	for (int i = 0; i < nval; ++i)
	{
		s += GetConcentration(*ppe[i], node_id);
	}
	s = s / nval;
	return s;
}

//-----------------------------------------------------------------------------
void FEPrescribedStressSensitiveConcentration::GetNodalValues(int nodelid, std::vector<double>& val)
{
	FENode* p_n = GetNodeSet()->Node(nodelid);
	int node_id = p_n->GetID() - 1;
	// get the stress
	mat3ds s = GetNodalStress(node_id);
	double c = GetNodalConcentration(node_id);
	double m_S = m_a0 + m_a / (1.0 + (exp(-(s.tr() - m_b) / stress0)));
	val[0] = max(c * m_S,0.0);
}