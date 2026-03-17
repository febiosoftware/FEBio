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
#include "FENewtonianFluid.h"
#include "FEFluid.h"
#include "FEViscConst.h"
#include <FECore/FEModel.h>

// define the material parameters
BEGIN_FECORE_CLASS(FENewtonianFluid, FEViscousFluid)
    ADD_PARAMETER(m_kappa, FE_RANGE_GREATER_OR_EQUAL(0.0), "kappa")->setUnits(UNIT_VISCOSITY)->setLongName("bulk viscosity");
    ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu"   )->setUnits(UNIT_VISCOSITY)->setLongName("shear viscosity");

// Optionally add strain-dependent normalized viscosity relations
    ADD_PROPERTY(m_kappahat, "kappahat", FEProperty::Optional)->SetLongName("normal. bulk viscosity vs strain");
    ADD_PROPERTY(m_muhat, "muhat", FEProperty::Optional)->SetLongName("normal. bulk viscosity vs strain");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FENewtonianFluid::FENewtonianFluid(FEModel* pfem) : FEViscousFluid(pfem)
{
    m_kappa = 0;
    m_mu = 0;
    m_kappahat = nullptr;
    m_muhat = nullptr;
}

//-----------------------------------------------------------------------------
//! initialization
bool FENewtonianFluid::Init()
{
    FEModel* pfem = GetFEModel();
    if (m_kappahat == nullptr) {
        m_kappahat = new FEViscConst(pfem);
    }
    m_kappahat->Init();
    if (m_muhat == nullptr) {
        m_muhat = new FEViscConst(pfem);
    }
    m_muhat->Init();
    return FEViscousFluid::Init();
}

//-----------------------------------------------------------------------------
void FENewtonianFluid::Serialize(DumpStream& ar)
{
    FEViscousFluid::Serialize(ar);
    
    if (ar.IsShallow()) return;
    
}

//-----------------------------------------------------------------------------
//! viscous stress
mat3ds FENewtonianFluid::Stress(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    
    mat3ds D = vt.RateOfDeformation();
    
    double mu = ShearViscosity(pt);
    double kappa = BulkViscosity(pt);
    
    mat3ds s = mat3dd(1.0)*(D.tr()*(kappa - 2.*mu/3.)) + D*(2*mu);
        
    return s;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to strain J
mat3ds FENewtonianFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& vt = *mp.ExtractData<FEFluidMaterialPoint>();
    
    mat3ds D = vt.RateOfDeformation();
    
    double dmudJ = m_mu*m_muhat->Tangent_NormalizedViscosity_Strain(mp);
    double dkappadJ = m_kappa*m_kappahat->Tangent_NormalizedViscosity_Strain(mp);
    
    mat3ds dsdJ = mat3dd(1.0)*(D.tr()*(dkappadJ - 2.*dmudJ/3.)) + D*(2*dmudJ);
        
    return dsdJ;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to rate of deformation tensor D
tens4ds FENewtonianFluid::Tangent_RateOfDeformation(FEMaterialPoint& mp)
{
    mat3dd I(1.0);
    double mu = ShearViscosity(mp);
    double kappa = BulkViscosity(mp);
    
    tens4ds c = dyad1s(I)*(kappa - 2.*mu/3.) + dyad4s(I)*(2*mu);
    
    return c;
}

//-----------------------------------------------------------------------------
//! dynamic shear viscosity
double FENewtonianFluid::ShearViscosity(FEMaterialPoint& mp)
{
    return m_mu*m_muhat->NormalizedViscosity(mp);
}

//! derivative of shear viscosity w.r.t. strain rate
double FENewtonianFluid::Tangent_ShearViscosity_StrainRate(FEMaterialPoint& mp)
{
	return 0.0;
}

//-----------------------------------------------------------------------------
//! bulk viscosity
double FENewtonianFluid::BulkViscosity(FEMaterialPoint& mp)
{
    return m_kappa*m_kappahat->NormalizedViscosity(mp);
}
