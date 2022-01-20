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
#include "FEVariableCriterion.h"
#include <FECore/FESolidDomain.h>
#include <FECore/FEMesh.h>

BEGIN_FECORE_CLASS(FEVariableCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_dof, "dof");
END_FECORE_CLASS();

FEVariableCriterion::FEVariableCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_dof = -1;
}

bool FEVariableCriterion::GetMaterialPointValue(FEMaterialPoint& mp, double& value)
{
	if (m_dof == -1) return false;

	FEElement* pe = mp.m_elem;
	if (pe == nullptr) return false;

	if ((mp.m_index < 0) || (mp.m_index >= pe->GaussPoints())) return false;

	FESolidDomain* dom = dynamic_cast<FESolidDomain*>(pe->GetMeshPartition());
	if (dom == nullptr) return false;

	FEMesh& mesh = *dom->GetMesh();

	double vn[FEElement::MAX_NODES];
	for (int i = 0; i < pe->Nodes(); ++i)
	{
		vn[i] = mesh.Node(pe->m_node[i]).get(m_dof);
	}

	value = pe->Evaluate(vn, mp.m_index);
	return true;
}
