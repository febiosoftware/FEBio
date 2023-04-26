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
}

//-----------------------------------------------------------------------------
mat3ds FEPrescribedStressSensitiveConcentration::GetStress(FEElement& m_elem, int node_id)
{
	mat3ds si[FEElement::MAX_INTPOINTS];
	mat3ds so[FEElement::MAX_NODES];

	for (int i = 0; i < m_elem.GaussPoints(); ++i) {
		FEMaterialPoint* pt = m_elem.GetMaterialPoint(i);
		FEElasticMaterialPoint* ep = pt->ExtractData<FEElasticMaterialPoint>();
		mat3dd I = mat3dd(1.0);
		if (ep)
			si[i] = ep->m_s + I * ep->m_p;
		else
			si[i].zero();
	}
	m_elem.project_to_nodes(si, so);
	return so[m_elem.FindNode(m_elem.m_node[node_id])];
}

//-----------------------------------------------------------------------------
double FEPrescribedStressSensitiveConcentration::GetEffectiveJacobian(FEElement& m_elem, int node_id)
{
	double ai[FEElement::MAX_INTPOINTS];
	double ao[FEElement::MAX_NODES];

	for (int i = 0; i < m_elem.GaussPoints(); ++i) {
		FEMaterialPoint* pt = m_elem.GetMaterialPoint(i);
		FEElasticMaterialPoint* ep = pt->ExtractData<FEElasticMaterialPoint>();
		FEBiphasicMaterialPoint* bp = pt->ExtractData<FEBiphasicMaterialPoint>();
		mat3dd I = mat3dd(1.0);
		if (ep && bp)
			ai[i] = ep->m_J - bp->m_phi0;
		else if (ep)
			ai[i] = ep->m_J;
		else
			ai[i] = 0.0;
	}
	m_elem.project_to_nodes(ai, ao);
	return ao[m_elem.FindNode(m_elem.m_node[node_id])];
}

//-----------------------------------------------------------------------------
double FEPrescribedStressSensitiveConcentration::GetConcentration(FEElement& m_elem, int node_id)
{
	double ai[FEElement::MAX_INTPOINTS];
	double ao[FEElement::MAX_NODES];

	for (int i = 0; i < m_elem.GaussPoints(); ++i) {
		FEMaterialPoint& pt = *m_elem.GetMaterialPoint(i);
		if (&pt)
			ai[i] = m_value(pt);
		else
			ai[i] = 0.0;
	}
	m_elem.project_to_nodes(ai, ao);
	return ao[m_elem.FindNode(m_elem.m_node[node_id])];
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
	double vol = 0.0;
	for (int i = 0; i < nval; ++i)
	{
		double v = GetMesh().ElementVolume(*ppe[i]);
		s += GetStress(*ppe[i], i) * v;
		vol += v;
	}
	s = s / (nval * vol);
	return s;
}

//-----------------------------------------------------------------------------
double FEPrescribedStressSensitiveConcentration::GetNodalEffectiveJacobian(int node_id)
{
	// get the mesh to which this surface belongs
	FENodeElemList& NEL = GetMesh().NodeElementList();
	int nval = NEL.Valence(node_id);
	FEElement** ppe = NEL.ElementList(node_id);

	// Get the elements that the node belongs to
	double J_eff = 0.0;
	double vol = 0.0;
	for (int i = 0; i < nval; ++i)
	{
		double v = GetMesh().ElementVolume(*ppe[i]);
		J_eff += GetEffectiveJacobian(*ppe[i], i) * v;
		vol += v;
	}
	J_eff = J_eff / vol;
	return J_eff;
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
	double vol = 0.0;
	for (int i = 0; i < nval; ++i)
	{
		double v = GetMesh().ElementVolume(*ppe[i]);
		s += GetConcentration(*ppe[i], i) * v;
		vol += v;
	}
	s = s / vol;
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
	
	double J_eff = GetNodalEffectiveJacobian(node_id);
	val[0] = c / J_eff * m_S;
}