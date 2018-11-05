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
    bool Init() override;

public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;
    
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;
    
	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt) override;

	//! Serialization
	void Serialize(DumpStream& ar) override;

protected:
	// integrated Fiber density
	void IntegrateFiberDensity();
    
protected:
    FEElasticFiberMaterial*     m_pFmat;    // pointer to fiber material
	FEFiberDensityDistribution* m_pFDD;     // pointer to fiber density distribution
	FEFiberIntegrationScheme*   m_pFint;    // pointer to fiber integration scheme

	double m_IFD;      // integrated fiber density

	DECLARE_FECORE_CLASS();
};
