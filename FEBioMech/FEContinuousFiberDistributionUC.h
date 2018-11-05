//
//  FEContinuousFiberDistributionUC.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEContinuousFiberDistributionUC__
#define __FEBioMech__FEContinuousFiberDistributionUC__

#include "FEUncoupledMaterial.h"
#include "FEElasticFiberMaterialUC.h"
#include "FEFiberDensityDistribution.h"
#include "FEFiberIntegrationScheme.h"
#include "FEFiberMaterialPoint.h"

//  This material is a container for a fiber material, a fiber density
//  distribution, and an integration scheme.
//
class FEContinuousFiberDistributionUC : public FEUncoupledMaterial
{
public:
    FEContinuousFiberDistributionUC(FEModel* pfem);
    ~FEContinuousFiberDistributionUC();
    
    // Initialization
    bool Init() override;
    
public:
	//! calculate stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;
    
	//! calculate tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;
    
	//! calculate deviatoric strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& pt) override;
    
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() override;

protected:
	// integrated Fiber density
	void IntegrateFiberDensity();
    
public:
    FEElasticFiberMaterialUC*   m_pFmat;    // pointer to fiber material
	FEFiberDensityDistribution* m_pFDD;     // pointer to fiber density distribution
	FEFiberIntegrationScheme*	m_pFint;    // pointer to fiber integration scheme
	double	m_IFD;	// integrated fiber distribution

	DECLARE_FECORE_CLASS();
};

#endif /* defined(__FEBioMech__FEContinuousFiberDistributionUC__) */
