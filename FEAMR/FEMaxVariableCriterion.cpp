/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FEMaxVariableCriterion.h"
#include <FECore/FESolidDomain.h>

BEGIN_FECORE_CLASS(FEMaxVariableCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_maxValue, "max_value");
	ADD_PARAMETER(m_dof, "dof");
END_FECORE_CLASS();

FEMaxVariableCriterion::FEMaxVariableCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_maxValue = 0.0;
	m_dof = -1;
}

bool FEMaxVariableCriterion::Check(FEElement& el, double& elemVal)
{
	if (m_dof == -1) return false;

	FESolidDomain* dom = dynamic_cast<FESolidDomain*>(el.GetMeshPartition());
	if (dom == nullptr) return false;

	FEMesh& mesh = *dom->GetMesh();
	double maxVal = -1e99;
	for (int i = 0; i < el.Nodes(); ++i)
	{
		double vi = mesh.Node(el.m_node[i]).get(m_dof);
		if (vi > maxVal) maxVal = vi;
	}
	elemVal = maxVal;
	return (elemVal >= m_maxValue);
}