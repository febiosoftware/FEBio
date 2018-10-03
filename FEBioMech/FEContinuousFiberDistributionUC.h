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
    bool Init();
    
public:
	//! calculate stress at material point
	mat3ds DevStress(FEMaterialPoint& pt);
    
	//! calculate tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt);
    
	//! calculate deviatoric strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& pt);
    
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();

protected:
	// integrated Fiber density
	void IntegrateFiberDensity();
    
public:
    FEPropertyT<FEElasticFiberMaterialUC>   m_pFmat;    // pointer to fiber material
	FEPropertyT<FEFiberDensityDistribution> m_pFDD;     // pointer to fiber density distribution
	FEPropertyT<FEFiberIntegrationScheme>	m_pFint;    // pointer to fiber integration scheme
	double	m_IFD;	// integrated fiber distribution

	DECLARE_FECORE_CLASS();
};

#endif /* defined(__FEBioMech__FEContinuousFiberDistributionUC__) */
