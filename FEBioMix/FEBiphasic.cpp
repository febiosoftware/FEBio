#include "FEBiphasic.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
// Material parameters for the FEBiphasic material
BEGIN_PARAMETER_LIST(FEBiphasic, FEMaterial)
	ADD_PARAMETER(m_phi0 , FE_RANGE_CLOSED(0.0, 1.0), "phi0");
	ADD_PARAMETER(m_rhoTw, FE_RANGE_GREATER_OR_EQUAL(0.0), "fluid_density");
    ADD_PARAMETER(m_tau  , FE_RANGE_GREATER_OR_EQUAL(0.0), "tau");
END_PARAMETER_LIST();

//============================================================================
// FEBiphasicMaterialPoint
//============================================================================
FEBiphasicMaterialPoint::FEBiphasicMaterialPoint(FEMaterialPoint* ppt) : FEMaterialPoint(ppt) {}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEBiphasicMaterialPoint::Copy()
{
	FEBiphasicMaterialPoint* pt = new FEBiphasicMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FEBiphasicMaterialPoint::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_p << m_gradp << m_gradpp << m_w << m_pa << m_phi0 << m_phi0p << m_phi0hat << m_Jp;
	}
	else
	{
		ar >> m_p >> m_gradp >> m_gradpp >> m_w >> m_pa >> m_phi0 >> m_phi0p >> m_phi0hat >> m_Jp;
	}

	FEMaterialPoint::Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEBiphasicMaterialPoint::Init()
{
	m_p = m_pa = 0;
	m_gradp = m_gradpp = vec3d(0,0,0);
	m_w = vec3d(0,0,0);
	m_phi0 = m_phi0p = 0;
	m_phi0hat = 0;
	m_Jp = 1;

	FEMaterialPoint::Init();
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
    m_tau = 0;

	// set material properties
	AddProperty(&m_pSolid, "solid"         );
	AddProperty(&m_pPerm , "permeability"  );
	AddProperty(&m_pSupp , "solvent_supply", 0);
    AddProperty(&m_pAmom , "active_supply" , 0);
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPoint* FEBiphasic::CreateMaterialPointData() 
{
	FEBiphasicMaterialPoint* pt = new FEBiphasicMaterialPoint(m_pSolid->CreateMaterialPointData());
	pt->m_phi0 = m_phi0;
	return pt;
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
