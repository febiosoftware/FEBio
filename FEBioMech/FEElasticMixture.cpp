/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FEElasticMixture.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
FEElasticMixtureMaterialPoint::FEElasticMixtureMaterialPoint() : FEMaterialPointArray(new FEElasticMaterialPoint)
{ 
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEElasticMixtureMaterialPoint::Copy()
{
	FEElasticMixtureMaterialPoint* pt = new FEElasticMixtureMaterialPoint;
	pt->m_w = m_w;
	pt->m_mp = m_mp;
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FEElasticMixtureMaterialPoint::Init()
{
	// allocate weight array
	m_w.resize(Components(), 1.0);

	// don't forget to initialize the base class
    FEMaterialPointArray::Init();
}

//-----------------------------------------------------------------------------
void FEElasticMixtureMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointArray::Serialize(ar);
	ar & m_w;
}

//=============================================================================
//								FEElasticMixture
//=============================================================================

BEGIN_FECORE_CLASS(FEElasticMixture, FEElasticMaterial)
	ADD_PROPERTY(m_pMat, "solid");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEElasticMixture::FEElasticMixture(FEModel* pfem) : FEElasticMaterial(pfem)
{
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEElasticMixture::CreateMaterialPointData() 
{ 
	FEElasticMixtureMaterialPoint* pt = new FEElasticMixtureMaterialPoint();
//	pt->SetName(m_pMat.GetName());
	int NMAT = Materials();
	for (int i=0; i<NMAT; ++i) 
	{
		FEMaterialPoint* pi = m_pMat[i]->CreateMaterialPointData();
		pt->AddMaterialPoint(pi);
	}
	return pt;
}

//-----------------------------------------------------------------------------
void FEElasticMixture::AddMaterial(FEElasticMaterial* pm) 
{ 
	m_pMat.push_back(pm); 
}

//-----------------------------------------------------------------------------
//! This function evaluates the stress at the material point by evaluating the
//! individual stress components. 
mat3ds FEElasticMixture::Stress(FEMaterialPoint& mp)
{
	FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
	vector<double>& w = pt.m_w;
	assert(w.size() == m_pMat.size());

	// get the elastic material point
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate stress
	mat3ds s; s.zero();
	for (int i=0; i < (int) m_pMat.size(); ++i)
	{
		FEMaterialPoint* mpi = pt.GetPointData(i);
		mpi->m_elem = mp.m_elem;
		mpi->m_index = mp.m_index;

		// copy the elastic material point data to the components
		FEElasticMaterialPoint& epi = *mpi->ExtractData<FEElasticMaterialPoint>();
		epi.m_rt = ep.m_rt;
		epi.m_r0 = mp.m_r0;// ep.m_r0;
		epi.m_F = ep.m_F;
		epi.m_J = ep.m_J;
        epi.m_v = ep.m_v;
        epi.m_a = ep.m_a;
        epi.m_L = ep.m_L;

		s += epi.m_s = m_pMat[i]->Stress(*mpi)*w[i];
	}

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEElasticMixture::Tangent(FEMaterialPoint& mp)
{
	FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
	vector<double>& w = pt.m_w;
	assert(w.size() == m_pMat.size());

	// get the elastic material point
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate elasticity tensor
	tens4ds c(0.);
	for (int i=0; i < (int) m_pMat.size(); ++i)
	{
		FEMaterialPoint* mpi = pt.GetPointData(i);
		mpi->m_elem = mp.m_elem;
		mpi->m_index = mp.m_index;

		// copy the elastic material point data to the components
		FEElasticMaterialPoint& epi = *mpi->ExtractData<FEElasticMaterialPoint>();
		epi.m_rt = ep.m_rt;
		epi.m_r0 = mp.m_r0;// ep.m_r0;
		epi.m_F = ep.m_F;
		epi.m_J = ep.m_J;
        epi.m_v = ep.m_v;
        epi.m_a = ep.m_a;
        epi.m_L = ep.m_L;

		c += m_pMat[i]->Tangent(*mpi)*w[i];
	}

	return c;
}

//-----------------------------------------------------------------------------
//! This function evaluates the stress at the material point by evaluating the
//! individual stress components.
double FEElasticMixture::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
	vector<double>& w = pt.m_w;
	assert(w.size() == m_pMat.size());
    
	// get the elastic material point
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// calculate strain energy density
	double sed = 0.0;
	for (int i=0; i < (int) m_pMat.size(); ++i)
	{
		FEMaterialPoint* mpi = pt.GetPointData(i);
		mpi->m_elem = mp.m_elem;
		mpi->m_index = mp.m_index;

		// copy the elastic material point data to the components
		FEElasticMaterialPoint& epi = *mpi->ExtractData<FEElasticMaterialPoint>();
		epi.m_rt = ep.m_rt;
		epi.m_r0 = mp.m_r0;// ep.m_r0;
		epi.m_F = ep.m_F;
		epi.m_J = ep.m_J;
        epi.m_v = ep.m_v;
        epi.m_a = ep.m_a;
        epi.m_L = ep.m_L;

		sed += m_pMat[i]->StrainEnergyDensity(*mpi)*w[i];
	}
    
	return sed;
}

//! specialized material points
void FEElasticMixture::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
    FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
    
    for (int i=0; i < (int) m_pMat.size(); ++i)
    {
        FEElasticMaterialPoint& epi = *pt.GetPointData(i)->ExtractData<FEElasticMaterialPoint>();
        epi.m_elem = mp.m_elem;
        FEMaterial* pmj = GetMaterial(i);
        pmj->UpdateSpecializedMaterialPoints(epi, tp);
    }
}

//! the density is the sum of the component densities
double FEElasticMixture::Density(FEMaterialPoint& mp)
{
	double d = 0.0;
	for (FEElasticMaterial* pe : m_pMat) d += pe->Density(mp);
	return d;
}
