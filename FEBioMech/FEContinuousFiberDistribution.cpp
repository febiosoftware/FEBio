/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FEContinuousFiberDistribution.h"

BEGIN_FECORE_CLASS(FEContinuousFiberDistribution, FEElasticMaterial)

	// material properties
	ADD_PROPERTY(m_pFmat, "fibers");
	ADD_PROPERTY(m_pFDD, "distribution");
	ADD_PROPERTY(m_pFint, "scheme");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEContinuousFiberDistribution::FEContinuousFiberDistribution(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_pFmat = 0;
	m_pFDD = 0;
	m_pFint = 0;
}

//-----------------------------------------------------------------------------
FEContinuousFiberDistribution::~FEContinuousFiberDistribution() {}

//-----------------------------------------------------------------------------
bool FEContinuousFiberDistribution::Init()
{
    // initialize base class
	if (FEElasticMaterial::Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! Serialization
void FEContinuousFiberDistribution::Serialize(DumpStream& ar)
{	
	FEElasticMaterial::Serialize(ar);
	if (ar.IsShallow()) return;
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FEContinuousFiberDistribution::Stress(FEMaterialPoint& mp)
{ 
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate stress
	mat3ds s; s.zero();

	// get the local coordinate systems
	mat3d Qt = GetLocalCS(mp).transpose();

	double IFD = 0.0;

	// obtain an integration point iterator
	FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&pt);
	if (it->IsValid())
	{
		do
		{
			// get the global fiber direction
			vec3d& N = it->m_fiber;

			// convert to local coordinates
			vec3d n0 = Qt*N;

			// rotate to local configuration to evaluate ellipsoidally distributed material coefficients
			double R = m_pFDD->FiberDensity(mp, n0);

			// integrate the fiber distribution
			IFD += R*it->m_weight;

			// calculate the stress
			double wn = it->m_weight;
			s += m_pFmat->Stress(pt, N)*(R*wn);
		}
		while (it->Next());
	}
	else IFD = 1.0;

	// don't forget to delete the iterator
	delete it;

	assert(IFD != 0.0);
	return s / IFD;
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FEContinuousFiberDistribution::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the local coordinate systems
	mat3d Qt = GetLocalCS(mp).transpose();

	// initialize stress tensor
	tens4ds c;
	c.zero();

	double IFD = 0.0;

	FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&pt);
	if (it->IsValid())
	{
		do
		{
			// get the global fiber direction
			vec3d& N = it->m_fiber;

			// convert to local
			vec3d n0 = Qt*N;

			// rotate to local configuration to evaluate ellipsoidally distributed material coefficients
			double R = m_pFDD->FiberDensity(mp, n0);

			// integrate the fiber distribution
			IFD += R*it->m_weight;

			// calculate the tangent
			c += m_pFmat->Tangent(mp, N)*(R*it->m_weight);
		}
		while (it->Next());
	}
	else IFD = 1.0;

	// don't forget to delete the iterator
	delete it;

	// we multiply by two to add contribution from other half-sphere
	assert(IFD != 0.0);
	return c / IFD;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
//double FEContinuousFiberDistribution::StrainEnergyDensity(FEMaterialPoint& pt) { return m_pFint->StrainEnergyDensity(pt); }
double FEContinuousFiberDistribution::StrainEnergyDensity(FEMaterialPoint& mp)
{ 
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the local coordinate systems
	mat3d Qt = GetLocalCS(mp).transpose();

	double IFD = 0.0;

	double sed = 0.0;
	FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&pt);
	if (it->IsValid())
	{
		do
		{
			// get fiber direction in global coordinate system
			vec3d& N = it->m_fiber;

			// convert to local coordinates
			vec3d n0 = Qt*N;

			// rotate to local configuration to evaluate ellipsoidally distributed material coefficients
			double R = m_pFDD->FiberDensity(mp, n0);

			// integrate the fiber distribution
			IFD += R*it->m_weight;

			// calculate the stress
			sed += m_pFmat->StrainEnergyDensity(mp, N)*(R*it->m_weight);
		}
		while (it->Next());
	}
	else IFD = 1.0;

	// don't forget to delete the iterator
	delete it;

	// we multiply by two to add contribution from other half-sphere
	assert(IFD != 0.0);
	return sed / IFD;
}
