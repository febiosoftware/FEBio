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
    bool Init();

public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
    
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);
    
	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt);

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();

	//! Serialization
	void Serialize(DumpStream& ar);

protected:
	// integrated Fiber density
	void IntegrateFiberDensity();
    
protected:
    FEPropertyT<FEElasticFiberMaterial>     m_pFmat;    // pointer to fiber material
	FEPropertyT<FEFiberDensityDistribution> m_pFDD;     // pointer to fiber density distribution
	FEPropertyT<FEFiberIntegrationScheme>   m_pFint;    // pointer to fiber integration scheme

	double m_IFD;      // integrated fiber density
};
