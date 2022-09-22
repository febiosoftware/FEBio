/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

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
FEMaterialPointData* FEContinuousFiberDistribution::CreateMaterialPointData()
{
	FEMaterialPointData* mp = FEElasticMaterial::CreateMaterialPointData();
	mp->SetNext(m_pFmat->CreateMaterialPointData());
    return mp;
}

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
    FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

	// calculate stress
	mat3ds s; s.zero();

	// get the local coordinate system
	mat3d Q = GetLocalCS(mp);
    
    double IFD = IntegratedFiberDensity(mp);

	// obtain an integration point iterator
	FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&mp);
	if (it->IsValid())
	{
		do
		{
			// get the fiber direction for that fiber distribution
			vec3d& N = it->m_fiber;

			// evaluate ellipsoidally distributed material coefficients
			double R = m_pFDD->FiberDensity(mp, N);
            
            // convert fiber to global coordinates
            vec3d n0 = Q*N;
            
			// calculate the stress
			double wn = it->m_weight;
			s += m_pFmat->FiberStress(mp, fp.FiberPreStretch(n0))*(R*wn);
		}
		while (it->Next());
	}

	// don't forget to delete the iterator
	delete it;

	// divide by IFD
	return s / IFD;
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FEContinuousFiberDistribution::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

	// get the local coordinate system
	mat3d Q = GetLocalCS(mp);
    
    double IFD = IntegratedFiberDensity(mp);

	// initialize stress tensor
	tens4ds c;
	c.zero();

	FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&mp);
	if (it->IsValid())
	{
		do
		{
            // get the fiber direction for that fiber distribution
            vec3d& N = it->m_fiber;
            
            // evaluate ellipsoidally distributed material coefficients
            double R = m_pFDD->FiberDensity(mp, N);
            
            // convert fiber to global coordinates
            vec3d n0 = Q*N;

			// calculate the tangent
			c += m_pFmat->FiberTangent(mp, fp.FiberPreStretch(n0))*(R*it->m_weight);
		}
		while (it->Next());
	}

	// don't forget to delete the iterator
	delete it;
    
	// divide by IFD
	return c / IFD;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
//double FEContinuousFiberDistribution::StrainEnergyDensity(FEMaterialPoint& pt) { return m_pFint->StrainEnergyDensity(pt); }
double FEContinuousFiberDistribution::StrainEnergyDensity(FEMaterialPoint& mp)
{ 
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

	// get the local coordinate system
	mat3d Q = GetLocalCS(mp);
    
    double IFD = IntegratedFiberDensity(mp);

	double sed = 0.0;
	FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&mp);
	if (it->IsValid())
	{
		do
		{
            // get the fiber direction for that fiber distribution
            vec3d& N = it->m_fiber;
            
            // evaluate ellipsoidally distributed material coefficients
            double R = m_pFDD->FiberDensity(mp, N);
            
            // convert fiber to global coordinates
            vec3d n0 = Q*N;

			// calculate the stress
			sed += m_pFmat->FiberStrainEnergyDensity(mp, fp.FiberPreStretch(n0))*(R*it->m_weight);
		}
		while (it->Next());
	}

	// don't forget to delete the iterator
	delete it;

	// divide by IFD
	return sed / IFD;
}

//-----------------------------------------------------------------------------
double FEContinuousFiberDistribution::IntegratedFiberDensity(FEMaterialPoint& mp)
{
	double IFD = 0;
	// NOTE: Pass nullptr to GetIterator to avoid issues with GK rule!
	FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(nullptr);
	if (it->IsValid())
	{
		do
		{
            // get the fiber direction for that fiber distribution
            vec3d& N = it->m_fiber;

			double R = m_pFDD->FiberDensity(mp, N);

			// integrate the fiber distribution
			IFD += R * it->m_weight;

		} while (it->Next());
	}

	// don't forget to delete the iterator
	delete it;

	// just in case
	if (IFD == 0.0) IFD = 1.0;

	return IFD;
}
