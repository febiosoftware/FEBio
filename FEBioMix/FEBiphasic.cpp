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
#include "FEBiphasic.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
// Material parameters for the FEBiphasic material
BEGIN_FECORE_CLASS(FEBiphasic, FEMaterial)
	ADD_PARAMETER(m_phi0 , FE_RANGE_CLOSED(0.0, 1.0), "phi0");
	ADD_PARAMETER(m_rhoTw, FE_RANGE_GREATER_OR_EQUAL(0.0), "fluid_density")->setUnits(UNIT_DENSITY);
    ADD_PARAMETER(m_tau  , FE_RANGE_GREATER_OR_EQUAL(0.0), "tau");

	// set material properties
	ADD_PROPERTY(m_pSolid, "solid", FEProperty::Required | FEProperty::TopLevel);
	ADD_PROPERTY(m_pPerm, "permeability");
	ADD_PROPERTY(m_pSupp, "solvent_supply", FEProperty::Optional);
	ADD_PROPERTY(m_pAmom, "active_supply", FEProperty::Optional);

    ADD_PROPERTY(m_Q, "mat_axis", FEProperty::Optional);

END_FECORE_CLASS();

//============================================================================
// FEBiphasicMaterialPoint
//============================================================================
FEBiphasicMaterialPoint::FEBiphasicMaterialPoint(FEMaterialPointData* ppt) : FEMaterialPointData(ppt) {}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEBiphasicMaterialPoint::Copy()
{
	FEBiphasicMaterialPoint* pt = new FEBiphasicMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FEBiphasicMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
	ar & m_p & m_gradp & m_gradpp;
	ar & m_w & m_pa & m_phi0 & m_phi0t & m_phi0p & m_phi0hat & m_Jp;
    ar & m_ss;
}

//-----------------------------------------------------------------------------
void FEBiphasicMaterialPoint::Init()
{
	m_p = m_pa = 0;
	m_gradp = m_gradpp = vec3d(0,0,0);
	m_w = vec3d(0,0,0);
	m_phi0 = m_phi0t = m_phi0p = 0;
	m_phi0hat = 0;
	m_Jp = 1;
    m_ss.zero();

	FEMaterialPointData::Init();
}

//============================================================================
// FEBiphasic
//============================================================================

//-----------------------------------------------------------------------------
//! FEBiphasic constructor

FEBiphasic::FEBiphasic(FEModel* pfem) : FEMaterial(pfem)
{ 
	m_rhoTw = 0; 
	m_phi0 = 0;
    m_tau = 0;

	m_pSolid = 0;
	m_pPerm = 0;
	m_pSupp = 0;
	m_pAmom = 0;
}

//-----------------------------------------------------------------------------
// initialize
bool FEBiphasic::Init()
{
    if (!m_pSolid->Init()) return false;
    if (!m_pPerm->Init()) return false;
    if (m_pSupp && !m_pSupp->Init()) return false;
    if (m_pAmom && !m_pAmom->Init()) return false;
    return FEMaterial::Init();
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPointData* FEBiphasic::CreateMaterialPointData()
{
	// create the solid material point
	FEMaterialPointData* ep = m_pSolid->CreateMaterialPointData();

    // create biphasic material point
    FEBiphasicMaterialPoint* pt = new FEBiphasicMaterialPoint(ep);
    
	return pt;
}

//-----------------------------------------------------------------------------
// update specialized material points
void FEBiphasic::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
    m_pSolid->UpdateSpecializedMaterialPoints(mp, tp);
    m_pPerm->UpdateSpecializedMaterialPoints(mp, tp);
    if (m_pSupp) m_pSupp->UpdateSpecializedMaterialPoints(mp, tp);
    if (m_pAmom) m_pAmom->UpdateSpecializedMaterialPoints(mp, tp);
}

//-----------------------------------------------------------------------------
//! Porosity in current configuration
double FEBiphasic::Porosity(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pet = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
	// relative volume
	double J = et.m_J;
	// porosity
//	double phiw = 1 - m_phi0/J;
	double phi0 = pet.m_phi0t;
	double phiw = 1 - phi0/J;
	// check for pore collapse
	// TODO: throw an error if pores collapse
	phiw = (phiw > 0) ? phiw : 0;
	
	return phiw;
}

//-----------------------------------------------------------------------------
//! The stress of a poro-elastic material is the sum of the fluid pressure
//! and the elastic stress. Note that this function is declared in the base class
//! so you do not have to reimplement it in a derived class, unless additional
//! pressure terms are required.

mat3ds FEBiphasic::Stress(FEMaterialPoint& mp)
{
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// calculate solid material stress
	mat3ds s = m_pSolid->Stress(mp);
	
	// add fluid pressure
	s.xx() -= pt.m_p;
	s.yy() -= pt.m_p;
	s.zz() -= pt.m_p;
	
	return s;
}

mat3ds FEBiphasic::SecantStress(FEMaterialPoint& mp)
{
    FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
    
    // calculate solid material stress
    mat3ds s = m_pSolid->SecantStress(mp);
    
    // add fluid pressure
    s.xx() -= pt.m_p;
    s.yy() -= pt.m_p;
    s.zz() -= pt.m_p;
    
    return s;
}

//-----------------------------------------------------------------------------
//! The tangent is the sum of the elastic tangent plus the fluid tangent. Note
//! that this function is declared in the base class, so you don't have to 
//! reimplement it unless additional tangent components are required.

tens4dmm FEBiphasic::Tangent(FEMaterialPoint& mp)
{
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// call solid tangent routine
	tens4dmm c = m_pSolid->Tangent(mp);
	
	// fluid pressure
	double p = pt.m_p;
	
	// adjust tangent for pressures
	double D[6][6] = {0};
	c.extract(D);
	
	D[0][0] -= -p;
	D[1][1] -= -p;
	D[2][2] -= -p;
	
	D[0][1] -= p; D[1][0] -= p;
	D[1][2] -= p; D[2][1] -= p;
	D[0][2] -= p; D[2][0] -= p;
	
	D[3][3] -= -p;
	D[4][4] -= -p;
	D[5][5] -= -p;
	
	return tens4dmm(D);
}

//-----------------------------------------------------------------------------
//! The tangent is the sum of the elastic tangent plus the fluid tangent. Note
//! that this function is declared in the base class, so you don't have to
//! reimplement it unless additional tangent components are required.

tens4dmm FEBiphasic::SecantTangent(FEMaterialPoint& mp)
{
    FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
    
    // call solid tangent routine
    tens4dmm c = m_pSolid->SecantTangent(mp);
    
    // fluid pressure
    double p = pt.m_p;
    
    // adjust tangent for pressures
    double D[6][6] = {0};
    c.extract(D);
    
    D[0][0] -= -p;
    D[1][1] -= -p;
    D[2][2] -= -p;
    
    D[0][1] -= p; D[1][0] -= p;
    D[1][2] -= p; D[2][1] -= p;
    D[0][2] -= p; D[2][0] -= p;
    
    D[3][3] -= -p;
    D[4][4] -= -p;
    D[5][5] -= -p;
    
    return tens4dmm(D);
}

//-----------------------------------------------------------------------------
//! actual fluid pressure (same as effective pressure here)

double FEBiphasic::Pressure(FEMaterialPoint& pt)
{
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
	return ppt.m_p;
}

//-----------------------------------------------------------------------------
//! Return the permeability tensor as a double array

void FEBiphasic::Permeability(double k[3][3], FEMaterialPoint& pt)

{
	mat3ds kt = m_pPerm->Permeability(pt);
	
	k[0][0] = kt.xx();
	k[1][1] = kt.yy();
	k[2][2] = kt.zz();
	k[0][1] = k[1][0] = kt.xy();
	k[1][2] = k[2][1] = kt.yz();
	k[2][0] = k[0][2] = kt.xz();
	
}

//-----------------------------------------------------------------------------
mat3ds FEBiphasic::Permeability(FEMaterialPoint& mp)
{
    return m_pPerm->Permeability(mp);
}

//-----------------------------------------------------------------------------
//! return tangent of permeability with strain
tens4dmm FEBiphasic::Tangent_Permeability_Strain(FEMaterialPoint& mp)
{
    return m_pPerm->Tangent_Permeability_Strain(mp);
}

//! return the material permeability property
mat3ds FEBiphasic::MaterialPermeability(FEMaterialPoint& mp, const mat3ds E)
{
    // Evaluate right Cauchy-Green tensor from E
    mat3ds C = mat3dd(1) + E*2;
    
    // Evaluate right stretch tensor U from C
    vec3d v[3];
    double lam[3];
    C.eigen2(lam, v);
    lam[0] = sqrt(lam[0]); lam[1] = sqrt(lam[1]); lam[2] = sqrt(lam[2]);
    mat3ds U = dyad(v[0])*lam[0] + dyad(v[1])*lam[1] + dyad(v[2])*lam[2];
    double J = lam[0]*lam[1]*lam[2];
    
    // temporarily replace F in material point with U
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3ds Ui = dyad(v[0])/lam[0] + dyad(v[1])/lam[1] + dyad(v[2])/lam[2];
    mat3d Fsafe = pt.m_F;
    double Jsafe = pt.m_J;
    pt.m_F = U;
    pt.m_J = J;
    
    // Evaluate hydraulic permeability
    mat3ds k = Permeability(mp);
    
    // Restore original F
    pt.m_F = Fsafe;
    pt.m_J = Jsafe;
    
    // Convert spatial permeability to material permeability
    mat3ds K = (Ui*k*Ui).sym()*J;

    return K;
}

//-----------------------------------------------------------------------------
//! calculate spatial tangent stiffness at material point, using secant method
tens4dmm FEBiphasic::SecantTangent_Permeability_Strain(FEMaterialPoint& mp)
{
    // extract the deformation gradient
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3d F = pt.m_F;
    double J = pt.m_J;
    mat3ds E = pt.Strain();
    mat3dd I(1);
    
    // calculate the material permeability at the current deformation gradient
    mat3ds K = MaterialPermeability(mp,E);
    
    // create deformation gradient increment
    double eps = 1e-9;
    vec3d e[3];
    e[0] = vec3d(1,0,0); e[1] = vec3d(0,1,0); e[2] = vec3d(0,0,1);
    tens4dmm Kmm;
    for (int k=0; k<3; ++k) {
        // evaluate incremental material permeability
        mat3ds dE = dyads(e[k], e[k])*(eps/2);
        mat3ds dK = (MaterialPermeability(mp,E+dE) - K)/eps;
        
        // evaluate the secant modulus
        Kmm(0,0,k,k) = dK.xx();
        Kmm(1,1,k,k) = dK.yy();
        Kmm(2,2,k,k) = dK.zz();
        Kmm(0,1,k,k) = Kmm(1,0,k,k) = dK.xy();
        Kmm(1,2,k,k) = Kmm(2,1,k,k) = dK.yz();
        Kmm(2,0,k,k) = Kmm(0,2,k,k) = dK.xz();
        for (int l=0; l<3; ++l) {
            if (l != k) {
                // evaluate incremental material permeability
                mat3ds dE = dyads(e[k], e[l])*(eps/2);
                mat3ds dK = (MaterialPermeability(mp,E+dE) - K)/eps;
                
                // evaluate the secant modulus
                Kmm(0,0,k,l) = Kmm(0,0,l,k) = dK.xx();
                Kmm(1,1,k,l) = Kmm(1,1,l,k) = dK.yy();
                Kmm(2,2,k,l) = Kmm(2,2,l,k) = dK.zz();
                Kmm(0,1,k,l) = Kmm(0,1,l,k) = Kmm(1,0,k,l) = Kmm(1,0,l,k) = dK.xy();
                Kmm(1,2,k,l) = Kmm(1,2,l,k) = Kmm(2,1,k,l) = Kmm(2,1,l,k) = dK.yz();
                Kmm(2,0,k,l) = Kmm(2,0,l,k) = Kmm(0,2,k,l) = Kmm(0,2,l,k) = dK.xz();
            }
        }
    }
    
    return Kmm.pp(F)/J;
}
