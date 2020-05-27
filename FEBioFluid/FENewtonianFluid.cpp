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
#include "FENewtonianFluid.h"
#include "FEFluid.h"
#include "FEBiphasicFSI.h"

// define the material parameters
BEGIN_FECORE_CLASS(FENewtonianFluid, FEViscousFluid)
    ADD_PARAMETER(m_kappa, FE_RANGE_GREATER_OR_EQUAL(0.0), "kappa");
    ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FENewtonianFluid::FENewtonianFluid(FEModel* pfem) : FEViscousFluid(pfem)
{
    m_kappa = 0;
    m_mu = 0;
}

//-----------------------------------------------------------------------------
//! viscous stress
mat3ds FENewtonianFluid::Stress(FEMaterialPoint& pt)
{
    FEMaterial* pm = GetFEModel()->GetMaterial(pt.m_elem->GetMatID());
    FEBiphasicFSI* bfsi = dynamic_cast<FEBiphasicFSI*> (pm);
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    FEBiphasicFSIMaterialPoint& bt = *pt.ExtractData<FEBiphasicFSIMaterialPoint>();
    
    mat3ds D = vt.RateOfDeformation();
    
    mat3ds s = mat3dd(1.0)*(D.tr()*(m_kappa - 2.*m_mu/3.)) + D*(2*m_mu);
    
    if (bfsi)
        s = s/bt.m_phif;
        
    return s;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to strain J
mat3ds FENewtonianFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    return mat3ds(0,0,0,0,0,0);
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to rate of deformation tensor D
tens4ds FENewtonianFluid::Tangent_RateOfDeformation(FEMaterialPoint& mp)
{
    FEMaterial* pm = GetFEModel()->GetMaterial(mp.m_elem->GetMatID());
    FEBiphasicFSI* bfsi = dynamic_cast<FEBiphasicFSI*> (pm);
    FEBiphasicFSIMaterialPoint& bt = *mp.ExtractData<FEBiphasicFSIMaterialPoint>();
    
    mat3dd I(1.0);
    tens4ds c = dyad1s(I)*(m_kappa - 2.*m_mu/3.) + dyad4s(I)*(2*m_mu);
    
    if (bfsi)
        c = c/bt.m_phif;
    
    return c;
}

//-----------------------------------------------------------------------------
//! dynamic shear viscosity
double FENewtonianFluid::ShearViscosity(FEMaterialPoint& mp)
{
    return m_mu;
}

//-----------------------------------------------------------------------------
//! bulke viscosity
double FENewtonianFluid::BulkViscosity(FEMaterialPoint& mp)
{
    return m_kappa;
}
