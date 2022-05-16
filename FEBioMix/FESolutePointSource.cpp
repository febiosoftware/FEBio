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
#include "FESolutePointSource.h"
#include "FESolute.h"
#include "FEMultiphasic.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FEElemElemList.h>
#include <algorithm>
#include <FEBioMech/FEElasticMaterialPoint.h>
#include <iostream>
#include <unordered_set>

BEGIN_FECORE_CLASS(FESolutePointSource, FEBodyLoad)
	ADD_PARAMETER(m_soluteId, "solute");
	ADD_PARAMETER(m_rate, "rate");
	ADD_PARAMETER(m_pos.x, "x");
	ADD_PARAMETER(m_pos.y, "y");
	ADD_PARAMETER(m_pos.z, "z");
END_FECORE_CLASS();

FESolutePointSource::FESolutePointSource(FEModel* fem) : FEBodyLoad(fem), m_search(&fem->GetMesh())
{
	m_dofC = -1;
	m_soluteId = -1;
	m_pos = vec3d(0, 0, 0);
	m_rate = 0.0;
}

vec3d FESolutePointSource::GetPosition() const
{
	return m_pos;
}

void FESolutePointSource::SetPosition(const vec3d& pos)
{
	m_pos = pos;
}

int FESolutePointSource::GetSoluteID() const
{
	return m_soluteId;
}

void FESolutePointSource::SetSoluteID(int soluteID)
{
	m_soluteId = soluteID;
}

double FESolutePointSource::GetRate() const
{
	return m_rate;
}

void FESolutePointSource::SetRate(double rate)
{
	m_rate = rate;
}

void FESolutePointSource::SetRadius(double radius)
{
	m_radius = radius;
}

void FESolutePointSource::SetAccumulateFlag(bool b) {
	m_accumulate = b;
}

void FESolutePointSource::SetAccumulateCAFlag(bool b) {
	m_accumulate_ca = b;
}


bool FESolutePointSource::Init()
{
	// see if the solute exists
	FEModel* fem = GetFEModel();

	bool bfound = false;
	int ndata = fem->GlobalDataItems();
	for (int i=0; i<ndata; ++i)
	{
		FESoluteData* soluteData = dynamic_cast<FESoluteData*>(fem->GetGlobalData(i));
		if (soluteData && (soluteData->GetID() == m_soluteId))
		{
			bfound = true;
			break;
		}
	}
	if (bfound == false) return false;

	// initialize octree search
	if (m_search.Init() == false) return false;

	// get the degree of freedom of the concentration
	m_dofC = fem->GetDOFIndex("concentration", m_soluteId - 1);

	return FEBodyLoad::Init();
}

// allow species to accumulate at the point source
void FESolutePointSource::Accumulate(double dc) {
	// find the element in which the point lies
	m_q[0] = m_q[1] = m_q[2] = 0.0;
	m_el = dynamic_cast<FESolidElement*>(m_search.FindElement(m_pos, m_q));
	if (m_el == nullptr) return;

	// make sure this element is part of a multiphasic domain
	FEDomain* dom = dynamic_cast<FEDomain*>(m_el->GetMeshPartition());
	FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
	if (mat == nullptr) return;

	const int nint = m_el->GaussPoints();

	// Make sure the material has the correct solute
	int solid = -1;
	int sols = mat->Solutes();
	for (int j = 0; j < sols; ++j)
	{
		int solj = mat->GetSolute(j)->GetSoluteID();
		if (solj == m_soluteId)
		{
			solid = j;
			break;
		}
	}
	if (solid == -1) return;

	m_rate = dc + m_rate;
	m_accumulate = true;
	m_accumulate_ca = true; //?? not needed?

}

void FESolutePointSource::Update()
{

}

//! Evaluate force vector
void FESolutePointSource::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	// get the domain in which this element resides
	FESolidDomain* dom = dynamic_cast<FESolidDomain*>(m_el->GetMeshPartition());
	FEMesh* mesh = dom->GetMesh();

	// get time increment
	double dt = tp.timeIncrement;

	// evaluate the shape functions at the position
	double H[FEElement::MAX_NODES];
	m_el->shape_fnc(H, m_q[0], m_q[1], m_q[2]);

	// evaluate the concentration at this point
	int neln = m_el->Nodes();
	double c[FEElement::MAX_NODES];
	for (int i = 0; i < neln; ++i) c[i] = mesh->Node(m_el->m_node[i]).get(m_dofC);
	double cx = m_el->evaluate(c, m_q[0], m_q[1], m_q[2]);
	double Ve = mesh->ElementVolume(*m_el);
	double v_rate = m_rate / Ve;

	// assemble the element load vector
	vector<double> fe(neln, 0.0);
	vector<int> lm(neln, -1);
	std::vector<FEMaterialPoint*> possible_ints;
	double total_elem = 0;
	FindIntInRadius(possible_ints, total_elem);
	if (possible_ints.size() == 0) {
		int nint = m_el->GaussPoints();
		double* w = m_el->GaussWeights();
		for (int n = 0; n < nint; ++n)
		{
			double* H_int = m_el->H(n);
			FEMaterialPoint* mp = m_el->GetMaterialPoint(n);
			double m_J = mp->m_J0;
			// loop over all nodes
			for (int j = 0; j < neln; ++j)
			{
				// get the value of the integrand for this node
				// -ve needed since external forces are subtracted from internal forces
				fe[j] += -v_rate * m_J * H[n] * H_int[j] * dt * w[n] * neln;
			}
		}

		//// get the LM vector
		for (int i = 0; i < neln; ++i) lm[i] = mesh->Node(m_el->m_node[i]).m_ID[m_dofC];
	}
	else {
		int nint = m_el->GaussPoints();
		double* w = m_el->GaussWeights();
		for (int n = 0; n < nint; ++n)
		{
			double* H_int = m_el->H(n);
			FEMaterialPoint* mp = m_el->GetMaterialPoint(n);
			double m_J = mp->m_J0;
			// loop over all nodes
			for (int j = 0; j < neln; ++j)
			{
				// get the value of the integrand for this node
				// -ve needed since external forces are subtracted from internal forces
				fe[j] += -v_rate * m_J * H[n] * H_int[j] * dt * w[n] * neln;
			}
		}

		////if (m_accumulate) { m_rate = 0; m_accumulate = false; }

		//// get the LM vector
		for (int i = 0; i < neln; ++i) lm[i] = mesh->Node(m_el->m_node[i]).m_ID[m_dofC];
	}
	// assemble into global vector

 	R.Assemble(lm, fe);
}

//! evaluate stiffness matrix
void FESolutePointSource::StiffnessMatrix(FELinearSystem& S, const FETimeInfo& tp)
{

	// get time increment
	double dt = tp.timeIncrement;

	// get the domain in which this element resides
	FESolidDomain* dom = dynamic_cast<FESolidDomain*>(m_el->GetMeshPartition());
	FEMesh* mesh = dom->GetMesh();

	// evaluate the shape functions at the position
	double H[FEElement::MAX_NODES];
	m_el->shape_fnc(H, m_q[0], m_q[1], m_q[2]);

	// evaluate the concentration at this point
	int neln = m_el->Nodes();
	int nint = m_el->GaussPoints();

	// assemble the element load vector
	FEElementMatrix ke(neln, neln); ke.zero();
	double* w = m_el->GaussWeights();
	double Ve = mesh->ElementVolume(*m_el);
	double v_rate = m_rate / Ve;
	for (int k = 0; k < nint; k++) {
		FEMaterialPoint* mp = m_el->GetMaterialPoint(k);
		double m_J = mp->m_J0;
		for (int i = 0; i < neln; ++i)
			for (int j = 0; j < neln; ++j)
			{
				ke[i][j] = H[i] * m_J * H[j] * v_rate * w[k] * dt;
			}
	}
	// get the LM vector
	vector<int> lm(neln, -1);
	for (int i = 0; i < neln; ++i) lm[i] = mesh->Node(m_el->m_node[i]).m_ID[m_dofC];

	// get the nodes
	vector<int> nodes(neln, -1);
	for (int i = 0; i < neln; ++i) nodes[i] = m_el->m_node[i];

	// assemble into global matrix
	ke.SetIndices(lm);
	ke.SetNodes(nodes);
	S.Assemble(ke);
}

void FESolutePointSource::FindIntInRadius(std::vector<FEMaterialPoint*>& possible_ints, double& total_elem) {

	// get element and set up buffers
	m_el = dynamic_cast<FESolidElement*>(m_search.FindElement(m_pos, m_q));
	std::unordered_set<FESolidElement*> visited;
	std::set<FESolidElement*> next;
	//std::vector<FEMaterialPoint*> possible_ints;
	visited.reserve(1000);
	possible_ints.reserve(500);

	//we will need to check the current element first
	next.insert(m_el);
	// create the element adjacency list.
	FEDomain* dom = dynamic_cast<FEDomain*>(m_el->GetMeshPartition());
	FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
	if (mat == nullptr) return;
	// calculate the element volume
	auto mesh = dom->GetMesh();
	// create the element-element list
	FEElemElemList EEL;
	EEL.Create(mesh);

	//while there are still elements to evaluate
	while (next.size()) {
		// get the element to be evaluated
		FESolidElement* cur = *next.begin();
		// remove the current element from the next buffer and add it to the visited buffer
		next.erase(next.begin());
		visited.insert(cur);
		// get the current element bounds
		std::vector<vec3d> cur_element_bounds;
		// add integration points within the radius
		bool int_flag = false;
		for (int i = 0; i < cur->GaussPoints(); i++)
		{
			FEMaterialPoint* mp = cur->GetMaterialPoint(i);
			vec3d disp = mp->m_r0 - m_pos;
			if (disp.norm() <= m_radius) {
				possible_ints.push_back(mp);
				int_flag = true;
			}
		}
		if (int_flag)
		{
			total_elem++;
		}

		// Add neighboring element to the next buffer as long as they haven't been visited.
		// get the global ID of the current element
		int cur_id = cur->GetID() - 1;
		// for each neighboring element
		for (int i = 0; i < EEL.NeighborSize(); i++)
		{
			if (EEL.Neighbor(cur_id, i))
			{
				// if that element has not been visited yet add it to the next list
				if (!visited.count(dynamic_cast<FESolidElement*>(EEL.Neighbor(cur_id, i))))
				{
					next.insert(dynamic_cast<FESolidElement*>(EEL.Neighbor(cur_id, i)));
				}
			}
		}
	}
}