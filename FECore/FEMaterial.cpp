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
	ADD_PROPERTY(m_pmap, "mat_axis", FEProperty::Optional);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMaterial::FEMaterial(FEModel* pfem) : FECoreBase(FEMATERIAL_ID), m_pfem(pfem)
{
	static int n = 1;
	m_pmap = nullptr;
}

//-----------------------------------------------------------------------------
FEMaterial::~FEMaterial()
{
	if (m_pmap) delete m_pmap; 
}

//-----------------------------------------------------------------------------
//! Get the model this material belongs to
FEModel* FEMaterial::GetFEModel()
{
	return m_pfem;
}

//-----------------------------------------------------------------------------
void FEMaterial::SetCoordinateSystemMap(FECoordSysMap* pmap)
{
	m_pmap = pmap;
}

//-----------------------------------------------------------------------------
FECoordSysMap* FEMaterial::GetCoordinateSystemMap()
{
	return m_pmap;
}

//-----------------------------------------------------------------------------
//! This function does nothing here. Derived classes will use this to set the 
//! local coordinate systems for material points.
void FEMaterial::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	
}

//-----------------------------------------------------------------------------
//! Initial material.
bool FEMaterial::Init()
{
	// initialize material axes
	if (m_pmap) m_pmap->Init();

	// initialize base class
	return FECoreBase::Init();
}

//-----------------------------------------------------------------------------
void FEMaterial::AddDomain(FEDomain* dom)
{
	m_domList.AddDomain(dom);
}

//-----------------------------------------------------------------------------
//! Store the material data to the archive
void FEMaterial::Serialize(DumpStream &ar)
{
	// We don't need to serialize material data for shallow copies.
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		// save the local coodinate system generator
		int nmap = (m_pmap ? 1 : 0);
		ar << nmap;
		if (m_pmap)
		{
			ar << m_pmap->GetTypeStr();
			m_pmap->Serialize(ar);
		}
	}
	else
	{
		// read the local cordinate system
		int nmap;
		ar >> nmap;
		if (m_pmap) delete m_pmap;
		m_pmap = 0;

		if (nmap)
		{
			FEModel& pfem = ar.GetFEModel();

			char sztype[64]={0};
			ar >> sztype;
			m_pmap = fecore_new<FECoordSysMap>(sztype, &pfem);
			m_pmap->Serialize(ar);
		}
	}

	// Save the material's parameters
	FECoreBase::Serialize(ar);
}
