#include "stdafx.h"
#include "FEMaterial.h"

REGISTER_SUPER_CLASS(FEMaterial, FEMATERIAL_ID);

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEMaterial, FECoreBase)
	ADD_PARAMETER(m_Q, "mat_axis");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMaterial::FEMaterial(FEModel* fem) : FECoreBase(fem, FEMATERIAL_ID)
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
