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
#include "FEFixedConcentration.h"
#include <FECore/FEModel.h>

BEGIN_FECORE_CLASS(FEFixedConcentration, FEBoundaryCondition)
	ADD_PARAMETER(m_sol, "solute", 0, "$(Solutes)");
	ADD_PROPERTY(m_nodeSet, "node_set", FEProperty::Reference);
END_FECORE_CLASS();

FEFixedConcentration::FEFixedConcentration(FEModel* fem) : FEFixedBC(fem)
{
	m_sol = 0;
}

bool FEFixedConcentration::Init()
{
	int nsol = m_sol - 1;
	if (nsol < 0) return false;

	FEModel* fem = GetFEModel();
	DOFS& dofs = fem->GetDOFS();
	m_dofs.clear();
	int nvar = dofs.GetVariableIndex("concentration");
	m_dofs.push_back(dofs.GetDOF(nvar, nsol));
	return FEFixedBC::Init();
}
