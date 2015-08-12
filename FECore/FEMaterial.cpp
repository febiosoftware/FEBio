// FEMaterial.cpp: implementation of the FEMaterial class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEMaterial.h"
#include <math.h>
#include <stdarg.h>
#include "FEModel.h"
#include "FECoreKernel.h"

//-----------------------------------------------------------------------------
MaterialError::MaterialError(const char* szfmt, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// print to file
	va_start(args, szfmt);
	vsprintf(m_szerr, szfmt, args);
	va_end(args);
}

//-----------------------------------------------------------------------------
FEMaterial::FEMaterial(FEModel* pfem) : FECoreBase(FEMATERIAL_ID), m_pfem(pfem)
{
	static int n = 1;
	m_szname[0] = 0;
	m_nID = -1;
	m_pmap = 0;
	m_pParent = 0;
	m_nRB = -1;
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
//! Sets the name of the material
void FEMaterial::SetName(const char* sz)
{ 
	strcpy(m_szname, sz); 
}

//-----------------------------------------------------------------------------
//! Return the name of the material
const char* FEMaterial::GetName()
{ 
	return m_szname; 
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
void FEMaterial::AddProperty(FEProperty* pp, const char* sz)
{
	pp->SetName(sz);
	m_Prop.push_back(pp);
}

//-----------------------------------------------------------------------------
int FEMaterial::Properties()
{
	int N = (int) m_Prop.size();
	int n = 0;
	for (int i=0; i<N; ++i) n += m_Prop[i]->size();
	return n;
}

//-----------------------------------------------------------------------------
FECoreBase* FEMaterial::GetProperty(int n)
{
	int N = (int) m_Prop.size();
	int m = 0;
	for (int i=0; i<N; ++i)
	{
		FEProperty* pm = m_Prop[i];
		int l = pm->size();
		if (m+l > n) return pm->get(n-m);
		m += l;
	}
	return 0;
}

//-----------------------------------------------------------------------------
int FEMaterial::FindPropertyIndex(const char* sz)
{
	int NP = (int) m_Prop.size();
	for (int i=0; i<NP; ++i)
	{
		const FEProperty* pm = m_Prop[i];
		if (pm && (strcmp(pm->GetName(), sz) == 0)) return i;
	}
	return -1;
}

//-----------------------------------------------------------------------------
bool FEMaterial::SetProperty(int i, FECoreBase* pb)
{
	FEProperty* pm = m_Prop[i];
	if (pm->IsType(pb)) { pm->SetProperty(pb); return true; }
	return false;
}

//-----------------------------------------------------------------------------
//! This function does nothing here. Derived classes will use this to set the 
//! local coordinate systems for material points.
void FEMaterial::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	
}

//-----------------------------------------------------------------------------
//! This function doesn't initialize anything but should be overridden in 
//! derived classes to initialize and check material parameters
void FEMaterial::Init()
{
	if (m_pmap) m_pmap->Init();
}

//-----------------------------------------------------------------------------
FEParam* FEMaterial::GetParameter(const ParamString& s)
{
	if (s.count() == 1) return FEMaterial::GetParameter(s);

	int NP = (int) m_Prop.size();
	for (int i=0; i<NP; ++i)
	{
		FEProperty* mp = m_Prop[i];
		if (s == mp->GetName()) return mp->GetParameter(s.next());
	}

	return 0;
}

//-----------------------------------------------------------------------------
//! Store the material data to the archive
void FEMaterial::Serialize(DumpFile &ar)
{
	// Save the material's parameters
	FEParamContainer::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_nID;
		ar << m_nRB;

		// save the local coodinate system generator
		int ntype = -1;
		if (m_pmap) ntype = m_pmap->m_ntype;
		else ntype = FE_MAP_NONE;
		assert(ntype != -1);
		ar << ntype;
		if (m_pmap) m_pmap->Serialize(ar);
	}
	else
	{
		ar >> m_nID;
		ar >> m_nRB;

		// read the local cordinate system
		int ntype;
		ar >> ntype;
		if (m_pmap) delete m_pmap;
		m_pmap = 0;
		assert(ntype != -1);
		FEModel* pfem = ar.GetFEModel();
		switch (ntype)
		{
		case FE_MAP_NONE    : m_pmap = 0; break;
		case FE_MAP_LOCAL   : m_pmap = new FELocalMap         (pfem); break;
		case FE_MAP_SPHERE  : m_pmap = new FESphericalMap     (pfem); break;
		case FE_MAP_CYLINDER: m_pmap = new FECylindricalMap   (pfem); break;
		case FE_MAP_VECTOR  : m_pmap = new FEVectorMap        (pfem); break;
		case FE_MAP_ANGLES  : m_pmap = new FESphericalAngleMap(pfem); break;
		}
		if (m_pmap) m_pmap->Serialize(ar);
	}

	// serialize all the material properties
	int NP = m_Prop.size();
	for (int i = 0; i<NP; ++i)
	{
		FEProperty* pmat = m_Prop[i];
		pmat->Serialize(ar);
	}
}
