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
#include "FEModelLoad.h"
#include "FESolver.h"

//-----------------------------------------------------------------------------
FEModelLoad::FEModelLoad(FEModel* pfem) : FEStepComponent(pfem), m_dof(pfem)
{
}

//-----------------------------------------------------------------------------
const FEDofList& FEModelLoad::GetDofList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
void FEModelLoad::Serialize(DumpStream& ar)
{
	FEStepComponent::Serialize(ar);
	ar & m_dof;
}

void FEModelLoad::PrepStep()
{

}

Matrix_Type FEModelLoad::PreferredMatrixType()
{
	FEParam* p = GetParameter("symmetric_stiffness");
	if (p && (p->type() == FE_PARAM_BOOL) && !p->value<bool>())
	{
		return REAL_UNSYMMETRIC;
	}
	return REAL_SYMMETRIC;
}

//-----------------------------------------------------------------------------
void FEModelLoad::LoadVector(FEGlobalVector& R)
{
	// base class does nothing
}

//-----------------------------------------------------------------------------
void FEModelLoad::StiffnessMatrix(FELinearSystem& LS)
{
	// base class does nothing.
}
