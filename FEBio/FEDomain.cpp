#include "stdafx.h"
#include "FEDomain.h"
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
// FEDomain
//-----------------------------------------------------------------------------
FEElement* FEDomain::FindElementFromID(int nid)
{
	for (int i=0; i<Elements(); ++i)
	{
		FEElement& el = ElementRef(i);
		if (el.m_nID == nid) return &el;
	}

	return 0;
}

//-----------------------------------------------------------------------------
void FEDomain::InitMaterialPointData()
{
	if (m_pMat == 0) return;

	for (int i=0; i<Elements(); ++i)
	{
		FEElement& el = ElementRef(i);
		for (int k=0; k<el.GaussPoints(); ++k) el.SetMaterialPointData(m_pMat->CreateMaterialPointData(), k);
	}
}
