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
#include "FEMaterial.h"
#include "DumpStream.h"

REGISTER_SUPER_CLASS(FEMaterial, FEMATERIAL_ID);

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEMaterial, FECoreBase)
	ADD_PARAMETER(m_Q, "mat_axis");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMaterial::FEMaterial(FEModel* fem) : FECoreBase(fem)
{
	static int n = 1;
	m_Q = mat3d::identity();
}

//-----------------------------------------------------------------------------
FEMaterial::~FEMaterial()
{
	for (size_t i = 0; i < m_param.size(); ++i) delete m_param[i];
	m_param.clear();
}

//-----------------------------------------------------------------------------
// evaluate local coordinate system at material point
mat3d FEMaterial::GetLocalCS(const FEMaterialPoint& mp)
{
	FEMaterial* parent = dynamic_cast<FEMaterial*>(GetParent());
	if (parent) {
		mat3d Qp = parent->GetLocalCS(mp); return Qp*m_Q(mp);
	}
	else return m_Q(mp);
}

//-----------------------------------------------------------------------------
//! Initial material.
bool FEMaterial::Init()
{
	// initialize base class
	return FECoreBase::Init();
}

//-----------------------------------------------------------------------------
void FEMaterial::AddDomain(FEDomain* dom)
{
	m_domList.AddDomain(dom);
}

//-----------------------------------------------------------------------------
FEDomainParameter* FEMaterial::FindDomainParameter(const std::string& paramName)
{
	for (int i = 0; i < m_param.size(); ++i)
	{
		FEDomainParameter* pi = m_param[i];
		if (pi->name() == paramName) return pi;
	}
	return nullptr;
}

//-----------------------------------------------------------------------------
void FEMaterial::AddDomainParameter(FEDomainParameter* p)
{
	assert(p);
	m_param.push_back(p);
}
