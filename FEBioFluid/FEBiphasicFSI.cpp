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
#include "FEBiphasicFSI.h"
#include <FECore/FECoreKernel.h>
#include <FECore/DumpStream.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEBiphasicFSI, FEMaterial)

ADD_PARAMETER(m_phi0 , FE_RANGE_CLOSED(0.0, 1.0), "phi0");

// material properties
ADD_PROPERTY(m_pSolid, "solid");
ADD_PROPERTY(m_pFluid, "fluid");
ADD_PROPERTY(m_pPerm, "permeability");

END_FECORE_CLASS();

//============================================================================
// FEFSIMaterialPoint
//============================================================================
FEBiphasicFSIMaterialPoint::FEBiphasicFSIMaterialPoint(FEMaterialPoint* pt) : FEMaterialPoint(pt) {}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEBiphasicFSIMaterialPoint::Copy()
{
    FEBiphasicFSIMaterialPoint* pt = new FEBiphasicFSIMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEBiphasicFSIMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPoint::Serialize(ar);
    ar & m_w & m_aw & m_Jdot & m_phis & m_phif & m_phi0 & m_gradphif & m_gradJ & m_Lw;
}

//-----------------------------------------------------------------------------
void FEBiphasicFSIMaterialPoint::Init()
{
    m_w = m_aw = m_gradphif = m_gradJ = vec3d(0,0,0);
    m_Jdot = m_phis = m_phif = m_phi0 = 0;
    m_Lw = mat3d(0.0);
    
    FEMaterialPoint::Init();
}

//============================================================================
// FEFluidFSI
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluidFSI constructor

FEBiphasicFSI::FEBiphasicFSI(FEModel* pfem) : FEMaterial(pfem)
{
    m_rhoTw = 0;
    m_phi0 = 0;
    
    m_pSolid = 0;
    m_pFluid = 0;
    m_pPerm = 0;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPoint* FEBiphasicFSI::CreateMaterialPointData()
{
    FEFluidMaterialPoint* fpt = new FEFluidMaterialPoint(m_pSolid->CreateMaterialPointData());
    FEBiphasicFSIMaterialPoint* bfpt = new FEBiphasicFSIMaterialPoint(fpt);
    
    bfpt->m_phi0 = m_phi0;
    
    return bfpt;
}

//-----------------------------------------------------------------------------
// initialize
bool FEBiphasicFSI::Init()
{
    m_rhoTw = m_pFluid->m_rhor;
    
    return FEMaterial::Init();
}

//-----------------------------------------------------------------------------
//! Porosity in current configuration
double FEBiphasicFSI::Porosity(FEMaterialPoint& pt)
{
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    
    // relative volume
    double J = et.m_J;
    double phiw = 1 - m_phi0/J;
    // check for pore collapse
    // TODO: throw an error if pores collapse
    // phiw cant be 0
    phiw = (phiw > 0) ? phiw : 1.0e-15;
    
    return phiw;
}

//-----------------------------------------------------------------------------
//! Solid Volume Frac (1 - porosity) in current configuration
double FEBiphasicFSI::SolidVolumeFrac(FEMaterialPoint& pt)
{
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    
    // relative volume
    double J = et.m_J;
    double phis = m_phi0/J;
    // check if phis is negative
    // TODO: throw an error if pores collapse
    phis = (phis >= 0) ? phis : 0;
    
    return phis;
}

//-----------------------------------------------------------------------------
//! The stress of a poro-elastic material is the sum of the fluid stress
//! and the elastic stress. Note that this function is declared in the base class
//! so you do not have to reimplement it in a derived class, unless additional
//! pressure terms are required.

mat3ds FEBiphasicFSI::Stress(FEMaterialPoint& mp)
{
    // calculate solid material stress
    mat3ds s = m_pSolid->Stress(mp);
    
    // add fluid stress
    s = s + m_pFluid->Stress(mp);
    
    return s;
}

//-----------------------------------------------------------------------------
//! The tangent is the sum of the elastic tangent plus the fluid tangent (it is 0). Note
//! that this function is declared in the base class, so you don't have to
//! reimplement it unless additional tangent components are required.

tens4ds FEBiphasicFSI::Tangent(FEMaterialPoint& mp)
{
    // call solid tangent routine
    tens4ds c = m_pSolid->Tangent(mp);
    return c;
}

//-----------------------------------------------------------------------------
//! Return the permeability tensor as a double array

void FEBiphasicFSI::Permeability(double k[3][3], FEMaterialPoint& pt)

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
mat3ds FEBiphasicFSI::Permeability(FEMaterialPoint& mp)
{
    return m_pPerm->Permeability(mp);
}

//-----------------------------------------------------------------------------
tens4dmm FEBiphasicFSI::Permeability_Tangent(FEMaterialPoint& mp)
{
    //Return 0 if permeability is 0 to avoid NAN
    if (m_pPerm->Permeability(mp).xx() == 0.0 && m_pPerm->Permeability(mp).xy() == 0.0 && m_pPerm->Permeability(mp).xz() == 0.0 && m_pPerm->Permeability(mp).yy() == 0.0 && m_pPerm->Permeability(mp).yz() == 0.0 && m_pPerm->Permeability(mp).zz() == 0.0)
        return tens4dmm(0.0);
    else
        return m_pPerm->Tangent_Permeability_Strain(mp);
}

//-----------------------------------------------------------------------------
mat3ds FEBiphasicFSI::InvPermeability(FEMaterialPoint& mp)
{
    //Return 0 for inverse permeability when permeability is set to 0.
    //Acts as if permeability if infinite.
    if (m_pPerm->Permeability(mp).xx() == 0.0 && m_pPerm->Permeability(mp).xy() == 0.0 && m_pPerm->Permeability(mp).xz() == 0.0 && m_pPerm->Permeability(mp).yy() == 0.0 && m_pPerm->Permeability(mp).yz() == 0.0 && m_pPerm->Permeability(mp).zz() == 0.0)
        return mat3ds(0.0);
    else
        return m_pPerm->Permeability(mp).inverse();
}

//-----------------------------------------------------------------------------
double FEBiphasicFSI::FluidDensity(FEMaterialPoint& mp)
{
    double rhoTf = TrueFluidDensity(mp);
    double phif = Porosity(mp);
    
    return rhoTf*phif;
}

//-----------------------------------------------------------------------------
double FEBiphasicFSI::SolidDensity(FEMaterialPoint& mp)
{
    double rhoTs = TrueSolidDensity(mp);
    double phis = SolidVolumeFrac(mp);
    
    return rhoTs*phis;
}
