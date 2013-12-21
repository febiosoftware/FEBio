#include "FEBiphasic.h"
#include "FECore/febio.h"

//-----------------------------------------------------------------------------
// Material parameters for the FEBiphasic material
BEGIN_PARAMETER_LIST(FEBiphasic, FEMaterial)
	ADD_PARAMETER(m_phi0, FE_PARAM_DOUBLE, "phi0");
	ADD_PARAMETER(m_rhoTw, FE_PARAM_DOUBLE, "fluid_density");
END_PARAMETER_LIST();

//============================================================================
// FEBiphasicMaterialPoint
//============================================================================
FEBiphasicMaterialPoint::FEBiphasicMaterialPoint(FEMaterialPoint* ppt) : FEMaterialPoint(ppt) {}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEBiphasicMaterialPoint::Copy()
{
	FEBiphasicMaterialPoint* pt = new FEBiphasicMaterialPoint(*this);
	if (m_pt) pt->m_pt = m_pt->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FEBiphasicMaterialPoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		dmp << m_p << m_gradp << m_w << m_pa << m_phi0 << m_phi0p << m_phi0hat;
	}
	else
	{
		dmp >> m_p >> m_gradp >> m_w >> m_pa >> m_phi0 >> m_phi0p >> m_phi0hat;
	}

	if (m_pt) m_pt->ShallowCopy(dmp, bsave);
}

//-----------------------------------------------------------------------------
void FEBiphasicMaterialPoint::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_p << m_gradp << m_w << m_pa << m_phi0 << m_phi0p << m_phi0hat;
	}
	else
	{
		ar >> m_p >> m_gradp >> m_w >> m_pa >> m_phi0 >> m_phi0p >> m_phi0hat;
	}

	if (m_pt) m_pt->Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEBiphasicMaterialPoint::Init(bool bflag)
{
	if (bflag)
	{
		m_p = m_pa = 0;
		m_gradp = vec3d(0,0,0);
		m_w = vec3d(0,0,0);
		m_phi0 = m_phi0p = 0;
		m_phi0hat = 0;
	}

	if (m_pt) m_pt->Init(bflag);
}

//============================================================================
// FEBiphasic
//============================================================================

//-----------------------------------------------------------------------------
//! FEBiphasic constructor

FEBiphasic::FEBiphasic(FEModel* pfem) : FEMaterial(pfem)
{ 
	m_rhoTw = 0; 
	m_phi0 = 0;

	m_pSolid = 0;
	m_pPerm = 0;
	m_pSupp = 0;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPoint* FEBiphasic::CreateMaterialPointData() 
{ 
	return new FEBiphasicMaterialPoint(m_pSolid->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
void FEBiphasic::Init()
{
	FEMaterial::Init();
	m_pSolid->SetParent(this); m_pSolid->Init();
	m_pPerm->SetParent(this); m_pPerm->Init();
	if (m_pSupp) { m_pSupp->SetParent(this); m_pSupp->Init(); }
	
	if (!INRANGE(m_phi0, 0.0, 1.0)) throw MaterialError("phi0 must be in the range 0 <= phi0 <= 1");
	if (m_rhoTw < 0) throw MaterialError("fluid_density must be positive");
}

//-----------------------------------------------------------------------------
//! A biphasic material has three properties
int FEBiphasic::Properties()
{
	return (m_pSupp ? 3 : 2);
}

//-----------------------------------------------------------------------------
//! return a pointer to a biphasic material property
FEMaterial* FEBiphasic::GetProperty(int i)
{
	switch (i)
	{
	case 0: return m_pSolid;
	case 1: return m_pPerm;
	case 2: return m_pSupp;
	}
	assert(false);
	return 0;
}

//-----------------------------------------------------------------------------
//! Find the index of a material property
int FEBiphasic::FindPropertyIndex(const char* szname)
{
	if (strcmp(szname, "solid"         ) == 0) return 0;
	if (strcmp(szname, "permeability"  ) == 0) return 1;
	if (strcmp(szname, "solvent_supply") == 0) return 2;
	return -1;
}

//-----------------------------------------------------------------------------
//! Set a material property
bool FEBiphasic::SetProperty(int n, FEMaterial* pm)
{
	switch(n)
	{
	case 0:
		{
			FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pm);
			if (pme) { m_pSolid = pme; return true; }
		}
		break;
	case 1: 
		{
			FEHydraulicPermeability* pmp = dynamic_cast<FEHydraulicPermeability*>(pm);
			if (pmp) { m_pPerm = pmp; return true; }
		}
		break;
	case 2:
		{
			FESolventSupply* pms = dynamic_cast<FESolventSupply*>(pm);
			if (pms) { m_pSupp = pms; return true; }
		}
		break;
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Porosity in current configuration
double FEBiphasic::Porosity(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pet = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
	// relative volume
	double J = et.m_J;
	// porosity
//	double phiw = 1 - m_phi0/J;
	double phi0 = pet.m_phi0;
	double phiw = 1 - phi0/J;
	// check for pore collapse
	// TODO: throw an error if pores collapse
	phiw = (phiw > 0) ? phiw : 0;
	
	return phiw;
}

//-----------------------------------------------------------------------------
//! The stress of a poro-elastic material is the sum of the fluid pressure
//! and the elastic stress. Note that this function is declared in the base class
//! so you do not have to reimplement it in a derived class, unless additional
//! pressure terms are required.

mat3ds FEBiphasic::Stress(FEMaterialPoint& mp)
{
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// calculate solid material stress
	mat3ds s = m_pSolid->Stress(mp);
	
	// add fluid pressure
	s.xx() -= pt.m_p;
	s.yy() -= pt.m_p;
	s.zz() -= pt.m_p;
	
	return s;
}

//-----------------------------------------------------------------------------
//! The tangent is the sum of the elastic tangent plus the fluid tangent. Note
//! that this function is declared in the base class, so you don't have to 
//! reimplement it unless additional tangent components are required.

tens4ds FEBiphasic::Tangent(FEMaterialPoint& mp)
{
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// call solid tangent routine
	tens4ds c = m_pSolid->Tangent(mp);
	
	// fluid pressure
	double p = pt.m_p;
	
	// adjust tangent for pressures
	double D[6][6] = {0};
	c.extract(D);
	
	D[0][0] -= -p;
	D[1][1] -= -p;
	D[2][2] -= -p;
	
	D[0][1] -= p; D[1][0] -= p;
	D[1][2] -= p; D[2][1] -= p;
	D[0][2] -= p; D[2][0] -= p;
	
	D[3][3] -= -p;
	D[4][4] -= -p;
	D[5][5] -= -p;
	
	return tens4ds(D);
}

//-----------------------------------------------------------------------------
//! Calculate fluid flux from the hydraulic permeability and the fluid pressure
//! gradient

vec3d FEBiphasic::Flux(FEMaterialPoint& pt)
{
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
	// pressure gradient
	vec3d gradp = ppt.m_gradp;
	
	// fluid flux w = -k*grad(p)
	mat3ds kt = m_pPerm->Permeability(pt);
	
	return -(kt*gradp);
}

//-----------------------------------------------------------------------------
//! actual fluid pressure (same as effective pressure here)

double FEBiphasic::Pressure(FEMaterialPoint& pt)
{
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
	return ppt.m_p;
}

//-----------------------------------------------------------------------------
//! Return the permeability tensor as a double array

void FEBiphasic::Permeability(double k[3][3], FEMaterialPoint& pt)

{
	mat3ds kt = m_pPerm->Permeability(pt);
	
	k[0][0] = kt.xx();
	k[1][1] = kt.yy();
	k[2][2] = kt.zz();
	k[0][1] = k[1][0] = kt.xy();
	k[1][2] = k[2][1] = kt.yz();
	k[2][0] = k[0][2] = kt.xz();
	
}

//-----------------------------------------------------------------------------
mat3ds FEBiphasic::Permeability(FEMaterialPoint& mp)
{
	return m_pPerm->Permeability(mp);
}

//-----------------------------------------------------------------------------
//! serialization
void FEBiphasic::Serialize(DumpFile &ar)
{
	// serialize material parameters
	FEMaterial::Serialize(ar);

	// serialize sub-materials
	int nSupp = 0;
	if (ar.IsSaving())
	{
		ar << m_pSolid->GetTypeStr();
		m_pSolid->Serialize(ar);

		ar << m_pPerm->GetTypeStr();
		m_pPerm->Serialize(ar);

		if (m_pSupp == 0) ar << nSupp;
		else
		{
			nSupp = 1;
			ar << nSupp;
			ar << m_pSupp->GetTypeStr();
			m_pSupp->Serialize(ar);
		}
	}
	else
	{
		char sz[256] = {0};

		ar >> sz;
		m_pSolid = dynamic_cast<FEElasticMaterial*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
		assert(m_pSolid);
		m_pSolid->Serialize(ar);
		m_pSolid->Init();

		ar >> sz;
		m_pPerm = dynamic_cast<FEHydraulicPermeability*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
		assert(m_pPerm);
		m_pPerm->Serialize(ar);
		m_pPerm->Init();

		ar >> nSupp;
		if (nSupp)
		{
			ar >> sz;
			m_pSupp = dynamic_cast<FESolventSupply*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
			assert(m_pSupp);
			m_pSupp->Serialize(ar);
			m_pSupp->Init();
		}
	}
}

//-----------------------------------------------------------------------------
FEParam* FEBiphasic::GetParameter(const ParamString& s)
{
	// see if this is a composite parameter
	if (s.count() == 1) return FEMaterial::GetParameter(s);

	// else find the component's parameter
	if      (s == "solid"       ) return m_pSolid->GetParameter(s.next());
	else if (s == "permeability") return m_pPerm ->GetParameter(s.next());
	else return 0;
}
