// FEMaterial.cpp: implementation of the FEMaterial class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEMaterial.h"
#include <math.h>
#include <stdarg.h>
#include "FECoreKernel.h"
#include "FEModelParam.h"

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
