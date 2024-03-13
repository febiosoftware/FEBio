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



#pragma once
#include "FEReactiveViscoelastic.h"
#include "FEReactiveFatigue.h"

//-----------------------------------------------------------------------------
// This material models damage in any reactive viscoelastic material.

class FERVEFatigueMaterial : public FEElasticMaterial
{
public:
    FERVEFatigueMaterial(FEModel* pfem);
    
public:
    //! calculate stress at material point
    mat3ds Stress(FEMaterialPoint& pt) override;
    
    //! calculate tangent stiffness at material point
    tens4ds Tangent(FEMaterialPoint& pt) override;
    
    //! calculate strain energy density at material point
    double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
    //! damage
    double Damage(FEMaterialPoint& pt);
        
    //! data initialization and checking
    bool Init() override;
    
    // returns a pointer to a new material point object
    FEMaterialPointData* CreateMaterialPointData() override;
    
    // get the elastic material
    FEElasticMaterial* GetElasticMaterial() override { return m_pRVE->GetElasticMaterial(); }
    
    //! specialized material points
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override
    {
        m_pRVE->UpdateSpecializedMaterialPoints(mp, tp);
    }
    
public:
    double StrongBondSED(FEMaterialPoint& pt) override;
    double WeakBondSED(FEMaterialPoint& pt) override;
    
public:
    FEReactiveViscoelasticMaterial* m_pRVE;     // reactive viscoelastic material

    DECLARE_FECORE_CLASS();
};
