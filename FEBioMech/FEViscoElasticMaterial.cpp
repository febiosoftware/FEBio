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
#include "FEViscoElasticMaterial.h"
#include "FEUncoupledMaterial.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/DumpStream.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEViscoElasticMaterial, FEElasticMaterial)

	// material parameters
	ADD_PARAMETER(m_t[0], "t1")->setUnits(UNIT_TIME);
	ADD_PARAMETER(m_t[1], "t2")->setUnits(UNIT_TIME);
	ADD_PARAMETER(m_t[2], "t3")->setUnits(UNIT_TIME);
	ADD_PARAMETER(m_t[3], "t4")->setUnits(UNIT_TIME);
	ADD_PARAMETER(m_t[4], "t5")->setUnits(UNIT_TIME);
	ADD_PARAMETER(m_t[5], "t6")->setUnits(UNIT_TIME);
	ADD_PARAMETER(m_g0  , "g0");
	ADD_PARAMETER(m_g[0], "g1");
	ADD_PARAMETER(m_g[1], "g2");
	ADD_PARAMETER(m_g[2], "g3");
	ADD_PARAMETER(m_g[3], "g4");
	ADD_PARAMETER(m_g[4], "g5");
	ADD_PARAMETER(m_g[5], "g6");

	// define the material properties
	ADD_PROPERTY(m_Base, "elastic");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEViscoElasticMaterialPoint::FEViscoElasticMaterialPoint(FEMaterialPointData* mp) : FEMaterialPointData(mp)
{
	m_sed = 0.0;
	m_sedp = 0.0;
}

//-----------------------------------------------------------------------------
//! Create a shallow copy of the material point data
FEMaterialPointData* FEViscoElasticMaterialPoint::Copy()
{
	FEViscoElasticMaterialPoint* pt = new FEViscoElasticMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! Initializes material point data.
void FEViscoElasticMaterialPoint::Init()
{
	// intialize data to zero
	m_Se.zero();
	m_Sep.zero();
	m_sed = 0.0;
    m_sedp = 0.0;
	for (int i=0; i<MAX_TERMS; ++i) {
		m_H[i].zero();
		m_Hp[i].zero();
        m_alpha[i] = m_alphap[i] = 1.0;
	}

    // don't forget to initialize the base class
    FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
//! Update material point data.
void FEViscoElasticMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	// the elastic stress stored in pt is the Cauchy stress.
	// however, we need to store the 2nd PK stress
	m_Sep = m_Se;
    m_sedp = m_sed;

	// copy previous data
	for (int i=0; i<MAX_TERMS; ++i) {
		m_Hp[i] = m_H[i];
        m_alphap[i] = m_alpha[i];
    }
    
    // don't forget to call the base class
    FEMaterialPointData::Update(timeInfo);
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FEViscoElasticMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);
	ar & m_Se;
	ar & m_Sep;
	ar & m_H & m_Hp;
    ar & m_sed & m_sedp;
    ar & m_alpha & m_alphap;
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

	m_Base = 0;
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
//! Create material point data for this material
FEMaterialPointData* FEViscoElasticMaterial::CreateMaterialPointData()
{
	return new FEViscoElasticMaterialPoint(m_Base->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FEViscoElasticMaterial::Stress(FEMaterialPoint& mp)
{
	double dt = GetFEModel()->GetTime().timeIncrement;
	if (dt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
    
	// get the elastic part
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the viscoelastic point data
	FEViscoElasticMaterialPoint& pt = *mp.ExtractData<FEViscoElasticMaterialPoint>();

	// Calculate the new elastic Cauchy stress
	mat3ds se = m_Base->Stress(mp);

	// pull-back to get PK2 stress
	mat3ds Se = pt.m_Se = ep.pull_back(se);

	// get elastic PK2 stress of previous timestep
	mat3ds Sep = pt.m_Sep;

	// calculate new history variables
	// terms are accumulated in S, the total PK2-stress
	mat3ds S = Se* m_g0(mp);
	double g, h;
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
	double dt = GetFEModel()->GetTime().timeIncrement;

	// calculate the spatial elastic tangent
	tens4ds C = m_Base->Tangent(pt);
	if (dt == 0.0) return C;

	// calculate the visco scale factor
	double f = m_g0(pt), g, h;
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
	// get the viscoelastic point data
	FEViscoElasticMaterialPoint& pt = *mp.ExtractData<FEViscoElasticMaterialPoint>();
    FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3d Fsafe = et.m_F; double Jsafe = et.m_J;

	// Calculate the new elastic strain energy density
	pt.m_sed = m_Base->StrainEnergyDensity(mp);
    double sed = pt.m_sed;
    
	double sedt = sed*m_g0(mp);
    if (SeriesStretchExponent(mp)) {
        // get the elastic point data and evaluate the right-stretch tensor
        for (int i=0; i<MAX_TERMS; ++i)
        {
            if (m_g[i] > 0) {
                mat3ds C = et.RightCauchyGreen();
                double l2[3], l[3];
                vec3d v[3];
                C.eigen2(l2, v);
                l[0] = sqrt(l2[0]); l[1] = sqrt(l2[1]); l[2] = sqrt(l2[2]);
                mat3ds Ua = dyad(v[0])*pow(l[0],pt.m_alpha[i])
                + dyad(v[1])*pow(l[1],pt.m_alpha[i]) + dyad(v[2])*pow(l[2],pt.m_alpha[i]);
                et.m_F = Ua; et.m_J = Ua.det();
                sedt += m_g[i]*m_Base->StrainEnergyDensity(mp);
            }
        }
    }
    else
        throw std::runtime_error("FEViscoElasticMaterial::strain energy density calculation did not converge!");
    
    et.m_F = Fsafe; et.m_J = Jsafe;

	// return the total strain energy density
	return sedt;
}

//-----------------------------------------------------------------------------
//! calculate exponent of right-stretch tensor in series spring
bool FEViscoElasticMaterial::SeriesStretchExponent(FEMaterialPoint& mp)
{
    const double errrel = 1e-6;
    const double almin = 0.001;
    const int maxiter = 50;
    // get the elastic point data and evaluate the right-stretch tensor
    FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();

    // get the right stretch tensor
    mat3ds C = et.RightCauchyGreen();
    double l2[3], l[3];
    vec3d v[3];
    C.eigen2(l2, v);
    l[0] = sqrt(l2[0]); l[1] = sqrt(l2[1]); l[2] = sqrt(l2[2]);
    mat3ds U = dyad(v[0])*l[0] + dyad(v[1])*l[1] + dyad(v[2])*l[2];
    double gamma = 0;
    for (int i=0; i<MAX_TERMS; ++i) gamma += m_g[i];
    
    // get the viscoelastic point data
    FEViscoElasticMaterialPoint& pt = *mp.ExtractData<FEViscoElasticMaterialPoint>();

    // use previous time solution as initial guess for the exponent
    mat3ds Se = pt.m_Se;
    mat3ds S = et.pull_back(et.m_s);
    double fmag = Se.dotdot(U);
    mat3d Fsafe = et.m_F; double Jsafe = et.m_J;
    for (int i=0; i<MAX_TERMS; ++i) {
        if (m_g[i] > 0) {
            double alpha = pt.m_alphap[i];
            bool done = false;
            int iter = 0;
            do {
                mat3ds Ua = dyad(v[0])*pow(l[0],alpha) + dyad(v[1])*pow(l[1],alpha) + dyad(v[2])*pow(l[2],alpha);
                et.m_F = Ua; et.m_J = Ua.det();
                mat3ds Sea = et.pull_back(m_Base->Stress(mp));
                double f = (Sea*m_g[i] - S + Se).dotdot(U);
                tens4ds Cea = et.pull_back(m_Base->Tangent(mp));
                mat3ds U2ap = dyad(v[0])*(pow(l[0],2*alpha)*log(l[0]))
                + dyad(v[1])*(pow(l[1],2*alpha)*log(l[1]))
                + dyad(v[2])*(pow(l[2],2*alpha)*log(l[2]));
                double fprime = (Cea.dot(U2ap)).dotdot(U)*m_g[i];
                if (fprime != 0) {
                    double dalpha = -f/fprime;
                    alpha += dalpha;
                    if (fabs(f) < errrel*fmag) done = true;
                    else if (fabs(dalpha) < errrel*fabs(alpha)) done = true;
                    else if (alpha > 1) { alpha = 1; done = true; }
                    else if (alpha < almin) { alpha = 0; done = true; }
                    else if (++iter > maxiter) done = true;
                }
                else
                    done = true;
            } while (!done);
            if (iter > maxiter) {
                et.m_F = Fsafe; et.m_J = Jsafe;
                return false;
            }
            pt.m_alpha[i] = alpha;
        }
    }
    et.m_F = Fsafe; et.m_J = Jsafe;
    return true;
}

