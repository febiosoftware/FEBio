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
#include "FEUncoupledElasticMixture.h"
#include "FECore/FECoreKernel.h"

// define the material parameters
 BEGIN_FECORE_CLASS(FEUncoupledElasticMixture, FEUncoupledMaterial)
	 ADD_PROPERTY(m_pMat, "solid");
 END_FECORE_CLASS();

//////////////////////////////////////////////////////////////////////
// Mixture of uncoupled elastic solids
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
FEUncoupledElasticMixture::FEUncoupledElasticMixture(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEUncoupledElasticMixture::CreateMaterialPointData() 
{ 
	FEElasticMixtureMaterialPoint* pt = new FEElasticMixtureMaterialPoint();
	int NMAT = Materials();
	for (int i=0; i<NMAT; ++i) pt->AddMaterialPoint(m_pMat[i]->CreateMaterialPointData());
	return pt;
}

//-----------------------------------------------------------------------------
bool FEUncoupledElasticMixture::Init()
{
	// NOTE: The calculation of K used to be the sum of all solid K's.
	//       But that doesn't follow the formulation and should be deprecated.
	//       Ideally, this function should be removed, but for backward compatiblity
	//       the old algorithm is retained (for now), only if the parent's K = 0. 
	//       Of course, if the user defined K for both the mixture and its components
	//       the behavior will be different. 
	if (m_K == 0.0)
	{
		for (int i=0; i < (int)m_pMat.size(); ++i) {
			m_pMat[i]->Init();
			m_K += m_pMat[i]->m_K;	// Sum up all the values of the bulk moduli
		}
	}

	return FEUncoupledMaterial::Init();
}

//-----------------------------------------------------------------------------
void FEUncoupledElasticMixture::AddMaterial(FEUncoupledMaterial* pm) 
{ 
	m_pMat.push_back(pm); 
}

//-----------------------------------------------------------------------------
mat3ds FEUncoupledElasticMixture::DevStress(FEMaterialPoint& mp)
{
	FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
	vector<double>& w = pt.m_w;
	assert(w.size() == m_pMat.size());

	// get the elastic material point
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate stress
	mat3ds s; s.zero();
	for (int i=0; i < (int)m_pMat.size(); ++i)
	{
		// copy the elastic material point data to the components
		FEElasticMaterialPoint& epi = *pt.GetPointData(i)->ExtractData<FEElasticMaterialPoint>();
		epi.m_rt = ep.m_rt;
		epi.m_r0 = ep.m_r0;
		epi.m_F = ep.m_F;
		epi.m_J = ep.m_J;
        epi.m_v = ep.m_v;
        epi.m_a = ep.m_a;
        epi.m_L = ep.m_L;

		s += epi.m_s = m_pMat[i]->DevStress(*pt.GetPointData(i))*w[i];
	}
    
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEUncoupledElasticMixture::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
	vector<double>& w = pt.m_w;
	assert(w.size() == m_pMat.size());

	// get the elastic material point
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate elasticity tensor
	tens4ds c(0.);
	for (int i=0; i < (int)m_pMat.size(); ++i)
	{
		// copy the elastic material point data to the components
		FEElasticMaterialPoint& epi = *pt.GetPointData(i)->ExtractData<FEElasticMaterialPoint>();
		epi.m_rt = ep.m_rt;
		epi.m_r0 = ep.m_r0;
		epi.m_F = ep.m_F;
		epi.m_J = ep.m_J;
        epi.m_v = ep.m_v;
        epi.m_a = ep.m_a;
        epi.m_L = ep.m_L;

		c += m_pMat[i]->DevTangent(*pt.GetPointData(i))*w[i];
	}
    
	return c;
}

//-----------------------------------------------------------------------------
double FEUncoupledElasticMixture::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
	vector<double>& w = pt.m_w;
	assert(w.size() == m_pMat.size());
    
	// get the elastic material point
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// calculate strain energy density
	double sed = 0.0;
	for (int i=0; i < (int)m_pMat.size(); ++i)
	{
		// copy the elastic material point data to the components
		FEElasticMaterialPoint& epi = *pt.GetPointData(i)->ExtractData<FEElasticMaterialPoint>();
		epi.m_rt = ep.m_rt;
		epi.m_r0 = ep.m_r0;
		epi.m_F = ep.m_F;
		epi.m_J = ep.m_J;
        epi.m_v = ep.m_v;
        epi.m_a = ep.m_a;
        epi.m_L = ep.m_L;

		sed += m_pMat[i]->DevStrainEnergyDensity(*pt.GetPointData(i))*w[i];
	}
    
	return sed;
}

//! the density is the sum of the component densities
double FEUncoupledElasticMixture::Density(FEMaterialPoint& mp)
{
	double d = 0.0;
	for (FEElasticMaterial* pe : m_pMat) d += pe->Density(mp);
	return d;
}
