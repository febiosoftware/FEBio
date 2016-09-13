//
//  FEContinuousFiberDistributionUC.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEContinuousFiberDistributionUC.h"

//-----------------------------------------------------------------------------
FEContinuousFiberDistributionUC::FEContinuousFiberDistributionUC(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
	m_IFD = 0.0;

	// set material properties
	AddProperty(&m_pFmat, "fibers"      );
	AddProperty(&m_pFDD , "distribution");
	AddProperty(&m_pFint, "scheme"      );
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
//! calculate stress at material point
mat3ds FEContinuousFiberDistributionUC::DevStress(FEMaterialPoint& mp)
{ 
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

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
			m_pFmat->SetFiberDirection(pt, n0);

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

			m_pFmat->SetFiberDirection(mp, n0e);

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
			m_pFmat->SetFiberDirection(mp, n0e);

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
