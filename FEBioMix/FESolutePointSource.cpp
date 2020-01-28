/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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

void FESolutePointSource::Update()
{
	// find the element in which the point lies
	m_q[0] = m_q[1] = m_q[2] = 0.0;
	m_el = dynamic_cast<FESolidElement*>(m_search.FindElement(m_pos, m_q));
	if (m_el == nullptr) return;

	// make sure this element is part of a multiphasic domain
	FEDomain* dom = dynamic_cast<FEDomain*>(m_el->GetMeshPartition());
	FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
	if (mat == nullptr) return;

	// Make sure the material has the correct solute
	int solid = -1;
	int sols = mat->Solutes();
	for (int j = 0; j<sols; ++j)
	{
		int sbmj = mat->GetSolute(j)->GetSoluteID();
		if (sbmj == m_soluteId)
		{
			solid = j;
			break;
		}
	}
	if (solid == -1) return;
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

	// assemble the element load vector
	vector<double> fe(neln, 0.0);
	for (int i = 0; i < neln; ++i)
	{
		fe[i] = -H[i] * m_rate * dt;
	}

	// get the LM vector
	vector<int> lm(neln, -1);
	for (int i = 0; i < neln; ++i) lm[i] = mesh->Node(m_el->m_node[i]).m_ID[m_dofC];

	// assemble into global vector
	R.Assemble(lm, fe);
}

//! evaluate stiffness matrix
void FESolutePointSource::StiffnessMatrix(FELinearSystem& S, const FETimeInfo& tp)
{
	return;

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

	// assemble the element load vector
	FEElementMatrix ke(neln, neln);
	for (int i = 0; i < neln; ++i)
		for (int j = 0; j < neln; ++j)
		{
			ke[i][j] = H[i] * H[j]*m_rate * dt;
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
