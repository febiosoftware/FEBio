//
//  FEContinuousFiberDistribution.h
//
//  Created by Gerard Ateshian on 11/17/13.
//

#pragma once
#include "FEElasticMaterial.h"
#include "FEElasticFiberMaterial.h"
#include "FEFiberDensityDistribution.h"
#include "FEFiberIntegrationScheme.h"
#include "FEFiberMaterialPoint.h"

//  This material is a container for a fiber material, a fiber density
//  distribution, and an integration scheme.
//
class FEContinuousFiberDistribution : public FEElasticMaterial
{
public:
    FEContinuousFiberDistribution(FEModel* pfem);
    ~FEContinuousFiberDistribution();
    
    // Initialization
    void Init();

public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) { return m_pFint->Stress(pt); }
    
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) { return m_pFint->Tangent(pt); }
    
	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt) { return m_pFint->StrainEnergyDensity(pt); }
    
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() {
        return new FEFiberMaterialPoint(m_pFint->CreateMaterialPointData());
    }
    
public:
    FEPropertyT<FEElasticFiberMaterial>     m_pFmat;    // pointer to fiber material
	FEPropertyT<FEFiberDensityDistribution> m_pFDD;     // pointer to fiber density distribution
	FEPropertyT<FEFiberIntegrationScheme>   m_pFint;    // pointer to fiber integration scheme
};
