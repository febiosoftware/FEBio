#include "stdafx.h"
#include "FEViscoElasticMaterial.h"
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEViscoElasticMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_t[0], FE_PARAM_DOUBLE, "t1");
	ADD_PARAMETER(m_t[1], FE_PARAM_DOUBLE, "t2");
	ADD_PARAMETER(m_t[2], FE_PARAM_DOUBLE, "t3");
	ADD_PARAMETER(m_t[3], FE_PARAM_DOUBLE, "t4");
	ADD_PARAMETER(m_t[4], FE_PARAM_DOUBLE, "t5");
	ADD_PARAMETER(m_t[5], FE_PARAM_DOUBLE, "t6");
	ADD_PARAMETER(m_g0  , FE_PARAM_DOUBLE, "g0");
	ADD_PARAMETER(m_g[0], FE_PARAM_DOUBLE, "g1");
	ADD_PARAMETER(m_g[1], FE_PARAM_DOUBLE, "g2");
	ADD_PARAMETER(m_g[2], FE_PARAM_DOUBLE, "g3");
	ADD_PARAMETER(m_g[3], FE_PARAM_DOUBLE, "g4");
	ADD_PARAMETER(m_g[4], FE_PARAM_DOUBLE, "g5");
	ADD_PARAMETER(m_g[5], FE_PARAM_DOUBLE, "g6");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Create a shallow copy of the material point data
FEMaterialPoint* FEViscoElasticMaterialPoint::Copy()
{
	FEViscoElasticMaterialPoint* pt = new FEViscoElasticMaterialPoint(*this);
	if (m_pt) pt->m_pt = m_pt->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! Initializes material point data.
void FEViscoElasticMaterialPoint::Init(bool bflag)
{
	FEElasticMaterialPoint& pt = *m_pt->ExtractData<FEElasticMaterialPoint>();
	if (bflag)
	{
		// intialize data to zero
		m_se.zero();
		m_Sep.zero();
		for (int i=0; i<MAX_TERMS; ++i) { m_H[i].zero(); m_Hp[i].zero(); };
	}
	else
	{
		// the elastic stress stored in pt is the Cauchy stress.
		// however, we need to store the 2nd PK stress
		m_Sep = pt.pull_back(m_se);

		// copy previous data
		for (int i=0; i<MAX_TERMS; ++i) m_Hp[i] = m_H[i];
	}

	// don't forget to intialize the nested data
	if (m_pt) m_pt->Init(bflag);
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FEViscoElasticMaterialPoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (m_pt) m_pt->ShallowCopy(dmp, bsave);

	if (bsave)
	{
		dmp << m_se;
		dmp << m_Sep;
		for (int i=0; i<MAX_TERMS; ++i) dmp << m_H[i] << m_Hp[i];
	}
	else
	{
		dmp >> m_se;
		dmp >> m_Sep;
		for (int i=0; i<MAX_TERMS; ++i) dmp >> m_H[i] >> m_Hp[i];
	}
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FEViscoElasticMaterialPoint::Serialize(DumpFile& ar)
{
	if (m_pt) m_pt->Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_se;
		ar << m_Sep;
		ar << (int) MAX_TERMS;
		for (int i=0; i<MAX_TERMS; ++i) ar << m_H[i] << m_Hp[i];
	}
	else
	{
		ar >> m_se;
		ar >> m_Sep;
		int n;
		ar >> n;
		assert(n == MAX_TERMS);
		for (int i=0; i<MAX_TERMS; ++i) ar >> m_H[i] >> m_Hp[i];
	}
}

//-----------------------------------------------------------------------------
//! constructor
FEViscoElasticMaterial::FEViscoElasticMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_g0 = 1;
	for (int i=0; i<MAX_TERMS; ++i)
	{
		m_t[i] = 1;
		m_g[i] = 0;
	}
	m_pBase = 0;
}

//-----------------------------------------------------------------------------
//! data initialization
//! \todo why does the base gets this material's parent?
void FEViscoElasticMaterial::Init()
{
	FEElasticMaterial::Init();
	if (m_pBase == 0) throw MaterialError("This material needs an elastic base.");
	m_pBase->SetParent(GetParent());
	m_pBase->Init();
}

//-----------------------------------------------------------------------------
//! This material only has one property
int FEViscoElasticMaterial::Properties()
{
	return 1;
}

//-----------------------------------------------------------------------------
FEMaterial* FEViscoElasticMaterial::GetProperty(int i)
{
	if (i == 0) return m_pBase;
	assert(false);
	return 0;
}

//-----------------------------------------------------------------------------
//! find a material property index ( returns <0 for error)
int FEViscoElasticMaterial::FindPropertyIndex(const char* szname)
{
	if (strcmp(szname, "elastic") == 0) return 0; else return -1;
}

//-----------------------------------------------------------------------------
//! set a material property (returns false on error)
bool FEViscoElasticMaterial::SetProperty(int i, FEMaterial* pm)
{
	if (i==0)
	{
		FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pm);
		if (pme && (dynamic_cast<FEUncoupledMaterial*>(pme) == 0))
		{ 
			SetBaseMaterial(pme);
			return true;
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPoint* FEViscoElasticMaterial::CreateMaterialPointData()
{ 
	return new FEViscoElasticMaterialPoint(m_pBase->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FEViscoElasticMaterial::Stress(FEMaterialPoint& mp)
{
	// get the elastic part
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the viscoelastic point data
	FEViscoElasticMaterialPoint& pt = *mp.ExtractData<FEViscoElasticMaterialPoint>();

	// Calculate the new elastic Cauchy stress
	pt.m_se = m_pBase->Stress(mp);

	// pull-back to get PK2 stress
	mat3ds Se = ep.pull_back(pt.m_se);

	// get elastic PK2 stress of previous timestep
	mat3ds Sep = pt.m_Sep;

	// calculate new history variables
	// terms are accumulated in S, the total PK2-stress
	mat3ds S = Se*m_g0;
	double dt = mp.dt, g, h;
	for (int i=0; i<MAX_TERMS; ++i)
	{
		g = exp(-dt/m_t[i]);
		h = (1 - g)/(dt/m_t[i]);

		pt.m_H[i] = pt.m_Hp[i]*g + (Se - Sep)*h;
		S += pt.m_H[i]*m_g[i];
	}

	// return the total Cauchy stress,
	// which is the push-forward of S
	return ep.push_forward(S);
}

//-----------------------------------------------------------------------------
//! Material tangent
tens4ds FEViscoElasticMaterial::Tangent(FEMaterialPoint& pt)
{
	// calculate the spatial elastic tangent
	tens4ds C = m_pBase->Tangent(pt);

	// calculate the visco scale factor
	double dt = pt.dt;
	double f = m_g0, g, h;
	for (int i=0; i<MAX_TERMS; ++i)
	{
		g = exp(-dt/m_t[i]);
		h = ( 1 - exp(-dt/m_t[i]) )/( dt/m_t[i] );
		f += m_g[i]*h; 
	}

	// multiply tangent with visco-factor
	return C*f;
}

//-----------------------------------------------------------------------------
//! Get a material parameter
FEParam* FEViscoElasticMaterial::GetParameter(const ParamString& s)
{
	if (s.count() == 1) return FEElasticMaterial::GetParameter(s);
	if (s == "elastic") return m_pBase->GetParameter(s.next());
	return 0;
}
