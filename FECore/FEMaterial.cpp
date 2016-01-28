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
bool MaterialError(const char* szfmt, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// create the message
	char szerr[512] = {0};
	va_start(args, szfmt);
	vsprintf(szerr, szfmt, args);
	va_end(args);

	return fecore_error(szerr);
}

//-----------------------------------------------------------------------------
FEProperty::FEProperty() : m_szname(nullptr), m_brequired(true) {}

//-----------------------------------------------------------------------------
FEProperty::~FEProperty(){}

//-----------------------------------------------------------------------------
//! Set the name of the property.
//! Note that the name is not copied so it must point to a static string.
FEProperty& FEProperty::SetName(const char* sz)
{ 
	m_szname = sz; 
	return *this; 
}

//-----------------------------------------------------------------------------
//! Return the name of this property
const char* FEProperty::GetName() const { return m_szname; }

//-----------------------------------------------------------------------------
void FEProperty::Write(DumpStream& ar, FEMaterial* pc)
{
	int nflag = (pc == 0 ? 0 : 1);
	ar << nflag;
	if (nflag)
	{
		ar << pc->GetTypeStr();
		pc->Serialize(ar);
	}
}

//-----------------------------------------------------------------------------
FEMaterial* FEProperty::Read(DumpStream& ar)
{
	int nflag = 0;
	FEMaterial* pm = 0;
	ar >> nflag;
	if (nflag)
	{
		char sz[256];
		ar >> sz;
		pm = fecore_new<FEMaterial>(FEMATERIAL_ID, sz, &ar.GetFEModel());
		pm->SetParent(GetParent());
		pm->Serialize(ar);

		// TODO: Do I really need to do this here?
		//pm->Init();
	}
	return pm;
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
void FEMaterial::AddProperty(FEProperty* pp, const char* sz, bool brequired)
{
	pp->SetName(sz);
	pp->m_brequired = brequired;
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
FEProperty* FEMaterial::FindProperty(const char* sz)
{
	int NP = (int) m_Prop.size();
	for (int i=0; i<NP; ++i)
	{
		FEProperty* pm = m_Prop[i];
		if (pm && (strcmp(pm->GetName(), sz) == 0)) return pm;
	}
	return nullptr;
}

//-----------------------------------------------------------------------------
bool FEMaterial::SetProperty(int i, FECoreBase* pb)
{
	FEProperty* pm = m_Prop[i];
	if (pm->IsType(pb))
	{ 
		pm->SetProperty(pb);
		FEMaterial* pmc = dynamic_cast<FEMaterial*>(pb);
		if (pmc) pmc->SetParent(this);
		return true; 
	}
	return false;
}

//-----------------------------------------------------------------------------
//! This function does nothing here. Derived classes will use this to set the 
//! local coordinate systems for material points.
void FEMaterial::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	
}

//-----------------------------------------------------------------------------
bool FEMaterial::Validate()
{
	// call base class first
	if (FEParamContainer::Validate() == false) return false;

	// check properties
	const int nprop = (int) m_Prop.size();
	for (int i=0; i<nprop; ++i) 
	{
		FEProperty* pi = m_Prop[i];
		if (pi)
		{
			if (pi->Validate() == false) return false;
		}
	}

	return true;

}

//-----------------------------------------------------------------------------
//! Initial material.
bool FEMaterial::Init()
{
	// check the parameter ranges
	if (Validate() == false) return false;

	// initialize material axes
	if (m_pmap) m_pmap->Init();

	// initialize all properties
	const int nprop = (int) m_Prop.size();
	for (int i=0; i<nprop; ++i) 
	{
		FEProperty* pi = m_Prop[i];
		if (pi)
		{
			if (pi->Init() == false) return false;
		}
		else return MaterialError("A nullptr was set for property i");
	}

	return true;
}

//-----------------------------------------------------------------------------
FEParam* FEMaterial::GetParameter(const ParamString& s)
{
	if (s.count() == 1) return FECoreBase::GetParameter(s);

	int NP = (int) m_Prop.size();
	for (int i=0; i<NP; ++i)
	{
		// get the property
		FEProperty* mp = m_Prop[i];

		// get the number of items in this property
		int nsize = mp->size();

		// If there is only one, we first compare the property title
		if ((nsize == 1) && (s == mp->GetName())) return mp->GetParameter(s.next());

		// for vector properties we try to match the "name" attribute
		for (int j=0; j<nsize; ++j)
		{
			FEMaterial* pm = dynamic_cast<FEMaterial*>(mp->get(j));
			if (pm && (s == pm->GetName())) return pm->GetParameter(s.next());
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
//! Store the material data to the archive
void FEMaterial::Serialize(DumpStream &ar)
{
	// We don't need to serialize material data for shallow copies.
	if (ar.IsShallow()) return;

	// Save the material's parameters
	FEParamContainer::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_nID;
		ar << m_nRB;

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
		ar >> m_nID;
		ar >> m_nRB;

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
			m_pmap = fecore_new<FECoordSysMap>(FECOORDSYSMAP_ID, sztype, &pfem);
			m_pmap->Serialize(ar);
		}
	}

	// serialize all the material properties
	int NP = m_Prop.size();
	for (int i = 0; i<NP; ++i)
	{
		FEProperty* pmat = m_Prop[i];
		pmat->SetParent(this);
		pmat->Serialize(ar);
	}
}

//-----------------------------------------------------------------------------
FEMaterial* FEMaterial::FindComponentByType(const char* sztype)
{
	// look at the properties first
	int NP = m_Prop.size();
	for (int i=0; i<NP; ++i)
	{
		// get the next property array
		FEProperty* pi = m_Prop[i];

		// loop over the properties in the array
		int n = pi->size();
		for (int j=0; j<n; ++j)
		{
			// get the next material
			FEMaterial* pmj = dynamic_cast<FEMaterial*>(pi->get(j));
			if (pmj)
			{
				// check the type string and return if match
				const char* sz = pmj->GetTypeStr();
				if (strcmp(sz, sztype)==0) return pmj;

				// If not matched, try recursively.
				FEMaterial* pmc = pmj->FindComponentByType(sztype);
				if (pmc) return pmc;
			}
		}
	}

	return 0;
}
