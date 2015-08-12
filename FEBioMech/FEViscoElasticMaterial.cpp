#include "stdafx.h"
#include "FEViscoElasticMaterial.h"
#include "FEUncoupledMaterial.h"
#include "FECore/FECoreKernel.h"

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
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! Initializes material point data.
void FEViscoElasticMaterialPoint::Init(bool bflag)
{
	FEElasticMaterialPoint& pt = *m_pNext->ExtractData<FEElasticMaterialPoint>();
	if (bflag)
	{
		// intialize data to zero
		m_se.zero();
		m_Sep.zero();
//        m_sed = 0.0;
//        m_sedp = 0.0;
		for (int i=0; i<MAX_TERMS; ++i) {
            m_H[i].zero();
            m_Hp[i].zero();
//            m_Hsed[i] = 0;
//            m_Hsedp[i] = 0;
        };
	}
	else
	{
		// the elastic stress stored in pt is the Cauchy stress.
		// however, we need to store the 2nd PK stress
		m_Sep = pt.pull_back(m_se);
//        m_sedp = m_sed;

		// copy previous data
		for (int i=0; i<MAX_TERMS; ++i) {
            m_Hp[i] = m_H[i];
//            m_Hsedp[i] = m_Hsed[i];
        }
	}

    // don't forget to initialize the base class
    FEMaterialPoint::Init(bflag);
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FEViscoElasticMaterialPoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (m_pNext) m_pNext->ShallowCopy(dmp, bsave);

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
	if (m_pNext) m_pNext->Serialize(ar);

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

	// define the material properties
	AddProperty(&m_Base, "elastic");
}

//-----------------------------------------------------------------------------
//! get the elastic base material \todo I want to call this GetElasticMaterial, but this name is being used
FEElasticMaterial* FEViscoElasticMaterial::GetBaseMaterial()
{ 
	return m_Base; 
}

//-----------------------------------------------------------------------------
//! Set the base material
void FEViscoElasticMaterial::SetBaseMaterial(FEElasticMaterial* pbase)
{ 
	m_Base = pbase; 
}

//-----------------------------------------------------------------------------
void FEViscoElasticMaterial::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	FEElasticMaterial::SetLocalCoordinateSystem(el, n, mp);
	FEElasticMaterial* pme2 = GetBaseMaterial();
	pme2->SetLocalCoordinateSystem(el, n, mp);
}

//-----------------------------------------------------------------------------
//! data initialization
//! \todo why does the base gets this material's parent?
void FEViscoElasticMaterial::Init()
{
	if (m_Base == 0) throw MaterialError("This material needs an elastic base.");
	FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPoint* FEViscoElasticMaterial::CreateMaterialPointData()
{ 
	return new FEViscoElasticMaterialPoint(m_Base->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FEViscoElasticMaterial::Stress(FEMaterialPoint& mp)
{
    if (mp.dt == 0) return mat3ds(0,0,0,0,0,0);
    
	// get the elastic part
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the viscoelastic point data
	FEViscoElasticMaterialPoint& pt = *mp.ExtractData<FEViscoElasticMaterialPoint>();

	// Calculate the new elastic Cauchy stress
	pt.m_se = m_Base->Stress(mp);

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
	tens4ds C = m_Base->Tangent(pt);
	if (pt.dt == 0.0) return C;

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
//! Strain energy density function
double FEViscoElasticMaterial::StrainEnergyDensity(FEMaterialPoint& mp)
{
/*    if (mp.dt == 0) return 0;
    
	// get the viscoelastic point data
	FEViscoElasticMaterialPoint& pt = *mp.ExtractData<FEViscoElasticMaterialPoint>();
    
	// Calculate the new elastic strain energy density
	pt.m_sed = m_pBase->StrainEnergyDensity(mp);
    double sed = pt.m_sed;
    
	// get elastic strain energy density of previous timestep
	double sedp = pt.m_sedp;
    
	// calculate new history variables
	// terms are accumulated in sedt, the total strain energy density
	double sedt = sed*m_g0;
	double dt = mp.dt, g, h;
	for (int i=0; i<MAX_TERMS; ++i)
	{
		g = exp(-dt/m_t[i]);
		h = (1 - g)/(dt/m_t[i]);
        
		pt.m_Hsed[i] = pt.m_Hsedp[i]*g + (sed - sedp)*h;
		sedt += pt.m_Hsed[i]*m_g[i];
	}
    
	// return the total strain energy density
	return sedt; */
    return 0;
}
