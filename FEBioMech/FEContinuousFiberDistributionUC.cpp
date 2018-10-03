#include "FEContinuousFiberDistributionUC.h"

BEGIN_FECORE_CLASS(FEContinuousFiberDistributionUC, FEUncoupledMaterial)
	// set material properties
	ADD_PROPERTY(m_pFmat, "fibers");
	ADD_PROPERTY(m_pFDD , "distribution");
	ADD_PROPERTY(m_pFint, "scheme");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEContinuousFiberDistributionUC::FEContinuousFiberDistributionUC(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
	m_IFD = 0.0;

	m_pFmat = 0;
	m_pFDD = 0;
	m_pFint = 0;
}

//-----------------------------------------------------------------------------
FEContinuousFiberDistributionUC::~FEContinuousFiberDistributionUC() {}

//-----------------------------------------------------------------------------
bool FEContinuousFiberDistributionUC::Init()
{
    m_K = m_pFmat->m_K;

	// initialize fiber integration scheme
	if (FEUncoupledMaterial::Init() == false) return false;

	// calculate the integrated fiber density
	IntegrateFiberDensity();

	return true;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPoint* FEContinuousFiberDistributionUC::CreateMaterialPointData() 
{
	return m_pFmat->CreateMaterialPointData();
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FEContinuousFiberDistributionUC::DevStress(FEMaterialPoint& mp)
{ 
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

	// calculate stress
	mat3ds s; s.zero();

	// get the element's local coordinate system
	mat3d QT = (pt.m_Q).transpose();

	// obtain an integration point iterator
	FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&pt);
	if (it->IsValid())
	{
		do
		{
			// set the fiber direction
			vec3d& n0 = it->m_fiber;
			fp.m_n0 = n0;

			// rotate to local configuration to evaluate ellipsoidally distributed material coefficients
			vec3d n0a = QT*n0;
			double R = m_pFDD->FiberDensity(n0a) / m_IFD;

			// calculate the stress
			double wn = it->m_weight;
			s += m_pFmat->DevStress(pt)*(R*wn);
		}
		while (it->Next());
	}

	// don't forget to delete the iterator
	delete it;

	return s;
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FEContinuousFiberDistributionUC::DevTangent(FEMaterialPoint& mp)
{ 
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

	// get the element's local coordinate system
	mat3d QT = (pt.m_Q).transpose();

	// initialize stress tensor
	tens4ds c;
	c.zero();

	FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&pt);
	if (it->IsValid())
	{
		do
		{
			// set fiber direction in global coordinate system
			vec3d& n0e = it->m_fiber;
			fp.m_n0 = n0e;

			// rotate to local configuration to evaluate ellipsoidally distributed material coefficients
			vec3d n0a = QT*n0e;
			double R = m_pFDD->FiberDensity(n0a) / m_IFD;

			// calculate the tangent
			c += m_pFmat->DevTangent(mp)*(R*it->m_weight);
		}
		while (it->Next());
	}

	// don't forget to delete the iterator
	delete it;

	// we multiply by two to add contribution from other half-sphere
	return c;
}

//-----------------------------------------------------------------------------
//! calculate deviatoric strain energy density
double FEContinuousFiberDistributionUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{ 
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

	// get the element's local coordinate system
	mat3d QT = (pt.m_Q).transpose();

	double sed = 0.0;
	FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&pt);
	if (it->IsValid())
	{
		do
		{
			// set fiber direction in global coordinate system
			vec3d& n0e = it->m_fiber;
			fp.m_n0 = n0e;

			// rotate to local configuration to evaluate ellipsoidally distributed material coefficients
			vec3d n0a = QT*n0e;
			double R = m_pFDD->FiberDensity(n0a) / m_IFD;

			// calculate the stress
			sed += m_pFmat->DevStrainEnergyDensity(mp)*(R*it->m_weight);
		}
		while (it->Next());
	}

	// don't forget to delete the iterator
	delete it;

	// we multiply by two to add contribution from other half-sphere
	return sed;
}

//-----------------------------------------------------------------------------
void FEContinuousFiberDistributionUC::IntegrateFiberDensity()
{
	m_IFD = 0;
	FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator();
	if (it->IsValid())
	{
		do
		{
			// set fiber direction in x-y plane of local coordinate system
			vec3d& n0a = it->m_fiber;

			// evaluate local fiber distribution
			double R = m_pFDD->FiberDensity(n0a);

			// integrate the fiber distribution
			m_IFD += R*it->m_weight;
		}
		while (it->Next());
	}

	// don't forget to delete the iterator
	delete it;
}
