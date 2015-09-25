#include "FEFluid.h"
#include "FECore/FECoreKernel.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEFluid, FEMaterial)
    ADD_PARAMETER2(m_rhor, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "density");
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
    m_vt = m_vp = m_at = vec3d(0,0,0);
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
void FEFluidMaterialPoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		dmp << m_p << m_L << m_J << m_Jp << m_gradJ << m_vt << m_vp << m_at << m_s;
	}
	else
	{
		dmp >> m_p >> m_L >> m_J >> m_Jp >> m_gradJ >> m_vt >> m_vp >> m_at >> m_s;
	}

	if (m_pNext) m_pNext->ShallowCopy(dmp, bsave);
}

//-----------------------------------------------------------------------------
void FEFluidMaterialPoint::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_p << m_L << m_J << m_Jp << m_gradJ << m_vt << m_vp << m_at << m_s;
	}
	else
	{
		ar >> m_p >> m_L >> m_J >> m_Jp >> m_gradJ >> m_vt >> m_vp >> m_at >> m_s;
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
