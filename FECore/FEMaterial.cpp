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
void FEProperty::Write(DumpFile& ar, FEMaterial* pc)
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
FEMaterial* FEProperty::Read(DumpFile& ar)
{
	int nflag = 0;
	FEMaterial* pm = 0;
	ar >> nflag;
	if (nflag)
	{
		char sz[256];
		ar >> sz;
		pm = fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel());
		pm->Serialize(ar);

		// TODO: Do I really need to do this here?
		pm->Init();
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
//! Initial material.
void FEMaterial::Init()
{
	// check the parameter ranges
	FEParameterList& pl = GetParameterList();
	const int nparam = pl.Parameters();
	list<FEParam>::iterator pi = pl.first();
	for (int i=0; i<nparam; ++i, ++pi)
	{
		FEParam& p = *pi;
		if (p.m_irange != FE_DONT_CARE)
		{
			if (p.is_inside_range() == false)
			{
				char szerr[256] = {0};
				const char* szname = p.m_szname;
				if (p.m_itype == FE_PARAM_INT)
				{
					switch (p.m_irange)
					{
					case FE_GREATER         : sprintf(szerr, "%s must be greater than %d"            , szname, p.m_imin); break;
					case FE_GREATER_OR_EQUAL: sprintf(szerr, "%s must be greater than or equal to %d", szname, p.m_imin); break;
					case FE_LESS            : sprintf(szerr, "%s must be less than %d"               , szname, p.m_imin); break;
					case FE_LESS_OR_EQUAL   : sprintf(szerr, "%s must be less than or equal to %d"   , szname, p.m_imin); break;
					case FE_OPEN            : sprintf(szerr, "%s must be in the open interval (%d, %d)"      , szname, p.m_imin, p.m_imax); break;
					case FE_CLOSED          : sprintf(szerr, "%s must be in the closed interval [%d, %d]"    , szname, p.m_imin, p.m_imax); break;
					case FE_LEFT_OPEN       : sprintf(szerr, "%s must be in the left-open interval (%d, %d]" , szname, p.m_imin, p.m_imax); break;
					case FE_RIGHT_OPEN      : sprintf(szerr, "%s must be in the right-open interval [%d, %d)", szname, p.m_imin, p.m_imax); break;
					case FE_NOT_EQUAL       : sprintf(szerr, "%s must not equal %d", szname, p.m_imin);
					default:
						sprintf(szerr, "%s has an invalid range");
					}
				}
				else if (p.m_itype == FE_PARAM_DOUBLE)
				{
					switch (p.m_irange)
					{
					case FE_GREATER         : sprintf(szerr, "%s must be greater than %lg"            , szname, p.m_dmin); break;
					case FE_GREATER_OR_EQUAL: sprintf(szerr, "%s must be greater than or equal to %lg", szname, p.m_dmin); break;
					case FE_LESS            : sprintf(szerr, "%s must be less than %lg"               , szname, p.m_dmin); break;
					case FE_LESS_OR_EQUAL   : sprintf(szerr, "%s must be less than or equal to %lg"   , szname, p.m_dmin); break;
					case FE_OPEN            : sprintf(szerr, "%s must be in the open interval (%lg, %lg)"      , szname, p.m_dmin, p.m_dmax); break;
					case FE_CLOSED          : sprintf(szerr, "%s must be in the closed interval [%lg, %lg]"    , szname, p.m_dmin, p.m_dmax); break;
					case FE_LEFT_OPEN       : sprintf(szerr, "%s must be in the left-open interval (%lg, %lg]" , szname, p.m_dmin, p.m_dmax); break;
					case FE_RIGHT_OPEN      : sprintf(szerr, "%s must be in the right-open interval [%lg, %lg)", szname, p.m_dmin, p.m_dmax); break;
					case FE_NOT_EQUAL       : sprintf(szerr, "%s must not equal %lg", szname, p.m_dmin);
					default:
						sprintf(szerr, "%s has an invalid range");
					}
				}
				else sprintf(szerr, "%s has an invalid range");

				// throw the error
				throw MaterialError(szerr);
			}
		}
	}

	// initialize material axes
	if (m_pmap) m_pmap->Init();

	// initialize all properties
	const int nprop = (int) m_Prop.size();
	for (int i=0; i<nprop; ++i) 
	{
		FEProperty* pi = m_Prop[i];
		if (pi)
		{
			if (pi->Init() == false)
			{
				// currently, the property will only return false if a required property was not defined
				throw MaterialError("This material requires the property %s", pi->GetName());
			}
		}
		else throw MaterialError("A nullptr was set for property i");
	}
}

//-----------------------------------------------------------------------------
FEParam* FEMaterial::GetParameter(const ParamString& s)
{
	if (s.count() == 1) return FECoreBase::GetParameter(s);

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
		case FE_MAP_POLAR   : m_pmap = new FEPolarMap         (pfem); break;
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
