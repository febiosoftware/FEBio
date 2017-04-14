#include "FEFluid.h"
#include "FECore/FECoreKernel.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEFluid, FEMaterial)
    ADD_PARAMETER2(m_rhor, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "density");
    ADD_PARAMETER2(m_k, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "k");
END_PARAMETER_LIST();

//============================================================================
// FEFluidMaterialPoint
//============================================================================
//-----------------------------------------------------------------------------
FEFluidMaterialPoint::FEFluidMaterialPoint()
{
    m_p = 0;
    m_L.zero();
    m_J = 1;
    m_Jdot = 0;
    m_vt = m_at = m_gradJ = vec3d(0,0,0);
    m_s.zero();
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEFluidMaterialPoint::Copy()
{
	FEFluidMaterialPoint* pt = new FEFluidMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FEFluidMaterialPoint::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_p << m_L << m_J << m_Jdot << m_gradJ << m_vt << m_at << m_s;
	}
	else
	{
		ar >> m_p >> m_L >> m_J >> m_Jdot >> m_gradJ >> m_vt >> m_at >> m_s;
	}

	if (m_pNext) m_pNext->Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEFluidMaterialPoint::Init()
{
	m_p = 0;
	m_L.zero();
	m_J = 1;
	m_Jdot = 0;
	m_vt = m_at = m_gradJ = vec3d(0,0,0);
	m_s.zero();
    
    // don't forget to initialize the base class
    FEMaterialPoint::Init();
}

//============================================================================
// FEFluid
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluid constructor

FEFluid::FEFluid(FEModel* pfem) : FEMaterial(pfem)
{ 
	m_rhor = 1;
    m_k = 1;

	// set material properties
	AddProperty(&m_pViscous ,"viscous"  );
}

//-----------------------------------------------------------------------------
//! returns a pointer to a new material point object
FEMaterialPoint* FEFluid::CreateMaterialPointData()
{
	FEFluidMaterialPoint* pt = new FEFluidMaterialPoint();
	return pt;
}

//-----------------------------------------------------------------------------
//! calculate current fluid density
double FEFluid::Density(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    return m_rhor/vt.m_J;
}

//-----------------------------------------------------------------------------
//! bulk modulus
double FEFluid::BulkModulus(FEMaterialPoint& mp)
{
    return m_k;
}

//-----------------------------------------------------------------------------
//! elastic pressure
double FEFluid::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double p = m_k*(1-fp.m_J);

    return p;
}

//-----------------------------------------------------------------------------
//! The stress of a fluid material is the sum of the fluid pressure
//! and the viscous stress.

mat3ds FEFluid::Stress(FEMaterialPoint& mp)
{
	// calculate solid material stress
	mat3ds s = m_pViscous->Stress(mp);
    
    double p = Pressure(mp);
	
	// add fluid pressure
	s.xx() -= p;
	s.yy() -= p;
	s.zz() -= p;
	
	return s;
}

//-----------------------------------------------------------------------------
//! The tangent of stress with respect to strain J of a fluid material is the
//! sum of the tangent of the fluid pressure and that of the viscous stress.

mat3ds FEFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    // get tangent of viscous stress
    mat3ds sJ = m_pViscous->Tangent_Strain(mp);
    
    // add tangent of fluid pressure
    sJ.xx() += m_k;
    sJ.yy() += m_k;
    sJ.zz() += m_k;
    
    return sJ;
}

//-----------------------------------------------------------------------------
//! calculate current fluid kinematic viscosity
double FEFluid::KinematicViscosity(FEMaterialPoint& mp)
{
    return m_pViscous->ShearViscosity(mp)/Density(mp);
}

//-----------------------------------------------------------------------------
//! calculate current acoustic speed
double FEFluid::AcousticSpeed(FEMaterialPoint& mp)
{
    double c = sqrt(BulkModulus(mp)/Density(mp));
    
    return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density (per referential volume)
double FEFluid::StrainEnergyDensity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double sed = m_k*pow(fp.m_J-1,2)/2;
    return sed;
}

//-----------------------------------------------------------------------------
//! calculate kinetic energy density (per referential volume)
double FEFluid::KineticEnergyDensity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double ked = m_rhor*(fp.m_vt*fp.m_vt)/2;
    return ked;
}

//-----------------------------------------------------------------------------
//! calculate energy density
double FEFluid::EnergyDensity(FEMaterialPoint& mp)
{
    return StrainEnergyDensity(mp) + KineticEnergyDensity(mp);
}
