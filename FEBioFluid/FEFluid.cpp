#include "FEFluid.h"
#include "FECore/FECoreKernel.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEFluid, FEMaterial)
    ADD_PARAMETER2(m_rhor, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "density");
    ADD_PARAMETER(m_bsupg, FE_PARAM_BOOL, "supg");
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
    m_Jp = 1;
    m_vt = m_vp = m_at = m_gradJ = vec3d(0,0,0);
    m_s.zero();
    m_L.zero();
    m_lapv = m_gdiv = vec3d(0,0,0);
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
		ar << m_p << m_L << m_J << m_Jp << m_gradJ << m_vt << m_vp << m_at << m_s << m_lapv << m_gdiv;
	}
	else
	{
		ar >> m_p >> m_L >> m_J >> m_Jp >> m_gradJ >> m_vt >> m_vp >> m_at >> m_s >> m_lapv >> m_gdiv;
	}

	if (m_pNext) m_pNext->Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEFluidMaterialPoint::Init(bool bflag)
{
	if (bflag)
	{
		m_p = 0;
        m_L.zero();
        m_J = 1;
        m_Jp = 1;
        m_vt = m_vp = m_at = m_gradJ = vec3d(0,0,0);
        m_s.zero();
        m_lapv = m_gdiv = vec3d(0,0,0);
	}

	if (m_pNext) m_pNext->Init(bflag);
    
    // don't forget to initialize the base class
    FEMaterialPoint::Init(bflag);
}

//============================================================================
// FEFluid
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluid constructor

FEFluid::FEFluid(FEModel* pfem) : FEMaterial(pfem)
{ 
	m_rhor = 1;
    m_bsupg = false;

	// set material properties
	AddProperty(&m_pElastic, "elastic"  );
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
//! The stress of a fluid material is the sum of the fluid pressure
//! and the viscous stress.

mat3ds FEFluid::Stress(FEMaterialPoint& mp)
{
	// calculate solid material stress
	mat3ds s = m_pViscous->Stress(mp);
    
    double p = m_pElastic->Pressure(mp);
	
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
    
    double dpdJ = m_pElastic->Tangent_Pressure_Strain(mp);
    
    // add tangent of fluid pressure
    sJ.xx() -= dpdJ;
    sJ.yy() -= dpdJ;
    sJ.zz() -= dpdJ;
    
    return sJ;
}

//-----------------------------------------------------------------------------
//! calculate current fluid kinematic viscosity
double FEFluid::KinematicViscosity(FEMaterialPoint& mp)
{
    return m_pViscous->DynamicViscosity(mp)/Density(mp);
}

//-----------------------------------------------------------------------------
//! calculate current acoustic speed
double FEFluid::AcousticSpeed(FEMaterialPoint& mp)
{
    double c = sqrt(m_pElastic->BulkModulus(mp)/Density(mp));
    
    return c;
}
