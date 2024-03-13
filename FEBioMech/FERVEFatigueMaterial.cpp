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
#include "FERVEFatigueMaterial.h"
#include "FEUncoupledMaterial.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERVEFatigueMaterial, FEElasticMaterial)
    // set material properties
    ADD_PROPERTY(m_pRVE , "viscoelastic");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FERVEFatigueMaterial::FERVEFatigueMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_pRVE  = nullptr;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FERVEFatigueMaterial::Init()
{
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>((FEElasticMaterial*)m_pRVE);
    if (m_pMat != nullptr)
    {
        feLogError("Viscoelastic material should not be of type uncoupled");
        return false;
    }
    
    FEReactiveFatigue* pmat = dynamic_cast<FEReactiveFatigue*>(m_pRVE->GetBaseMaterial());
    if (pmat == nullptr)
    {
        feLogError("Elastic material should be of type Reactive fatigue");
        return false;
    }
    
    return m_pRVE->Init();
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FERVEFatigueMaterial::CreateMaterialPointData()
{
    FEReactiveViscoelasticMaterialPoint* pt = new FEReactiveViscoelasticMaterialPoint();
    // create damage materal point for strong bond (base) material
    FEReactiveFatigueMaterialPoint* pbase = new FEReactiveFatigueMaterialPoint(m_pRVE->GetBaseMaterial()->CreateMaterialPointData());
    pt->AddMaterialPoint(new FEMaterialPoint(pbase));
    
    // create materal point for weak bond material
    FEReactiveVEMaterialPoint* pbond = new FEReactiveVEMaterialPoint(m_pRVE->GetBondMaterial()->CreateMaterialPointData());
    pt->AddMaterialPoint(new FEMaterialPoint(pbond));
    
    return pt;
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FERVEFatigueMaterial::Stress(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the stress
    mat3ds s = m_pRVE->Stress(pt);
    
    // return damaged stress
    return s*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FERVEFatigueMaterial::Tangent(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the tangent
    tens4ds c = m_pRVE->Tangent(pt);
    
    // return damaged tangent
    return c*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FERVEFatigueMaterial::StrainEnergyDensity(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the strain energy density
    double sed = m_pRVE->StrainEnergyDensity(pt);
    
    // return damaged sed
    return sed*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate strong bond strain energy density at material point
double FERVEFatigueMaterial::StrongBondSED(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the strain energy density
    double sed = m_pRVE->StrongBondSED(pt);
    
    // return damaged sed
    return sed*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate weak bond strain energy density at material point
double FERVEFatigueMaterial::WeakBondSED(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the strain energy density
    double sed = m_pRVE->WeakBondSED(pt);
    
    // return damaged sed
    return sed*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate damage at material point
double FERVEFatigueMaterial::Damage(FEMaterialPoint& pt)
{
    // get the reactive viscoelastic base material point
    FEMaterialPoint* pb = m_pRVE->GetBaseMaterialPoint(pt);
    
    // get the reactive fatigue material point data
    FEReactiveFatigueMaterialPoint& pd = *pb->ExtractData<FEReactiveFatigueMaterialPoint>();

    return pd.m_D;
}
