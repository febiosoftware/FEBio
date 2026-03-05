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
#include "FENewtonianThermoFluid.h"
#include <FEBioFluid/FEFluid.h>
#include <FEBioFluid/FEBiphasicFSI.h>
#include "FEThermoFluid.h"
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FENewtonianThermoFluid, FEThermoViscousFluid)
// properties
    ADD_PROPERTY(m_prop[0], "khat" ,FEProperty::Optional)->SetLongName("normalized bulk viscosity");
    ADD_PROPERTY(m_prop[1], "muhat",FEProperty::Optional)->SetLongName("normalized shear viscosity");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FENewtonianThermoFluid::FENewtonianThermoFluid(FEModel* pfem) : FEThermoViscousFluid(pfem)
{
}

//-----------------------------------------------------------------------------
//! initialization
bool FENewtonianThermoFluid::Init()
{
    if (m_prop[0]) m_prop[0]->Init();
    if (m_prop[1]) m_prop[1]->Init();
    
    // store an invariant copy of the Newtonian fluid viscosities
    m_kappa = m_NF->m_kappa;
    m_mu = m_NF->m_mu;
    
    return FEViscousFluid::Init();
}

//-----------------------------------------------------------------------------
//! viscous stress
mat3ds FENewtonianThermoFluid::Stress(FEMaterialPoint& pt)
{
    // evaluate thermal dependence of viscosities
    double kappa = m_kappa*ThermoProp(pt,0);
    double mu = m_mu*ThermoProp(pt,1);
    
    // set these viscosities in the Newtonian fluid object
    m_NF->m_kappa = kappa;
    m_NF->m_mu = mu;
    
    // evaluate the stress in the Newtonian fluid using these viscosities
    return m_NF->Stress(pt);
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to strain J
mat3ds FENewtonianThermoFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    return m_NF->Tangent_Strain(mp);
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to rate of deformation tensor D
tens4ds FENewtonianThermoFluid::Tangent_RateOfDeformation(FEMaterialPoint& mp)
{
    // evaluate thermal dependence of viscosities
    double kappa = m_kappa*ThermoProp(mp,0);
    double mu = m_mu*ThermoProp(mp,1);
    
    // set these viscosities in the Newtonian fluid object
    m_NF->m_kappa = kappa;
    m_NF->m_mu = mu;
    
    return m_NF->Tangent_RateOfDeformation(mp);
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to temperature T
mat3ds FENewtonianThermoFluid::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& vt = *mp.ExtractData<FEFluidMaterialPoint>();
    
    mat3ds D = vt.RateOfDeformation();
    
    double dkappa = m_kappa*TangentThermoProp(mp,0);
    double dmu = m_mu*TangentThermoProp(mp,1);
    
    mat3ds ds = mat3dd(1.0)*(D.tr()*(dkappa - 2.*dmu/3.)) + D*(2*dmu);
        
    return ds;
}
