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
#include "FEFiberIntegrationSchemeUC.h"
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
    void Init();
     
public:
	//! get the number of material properties
	int MaterialProperties();
    
	//! get a specific material property
	FEProperty* GetMaterialProperty(int i);
    
public:
	//! calculate stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) { return m_pFint->DevStress(pt); }
    
	//! calculate tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) { return m_pFint->DevTangent(pt); }
    
	//! calculate deviatoric strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& pt) { return m_pFint->DevStrainEnergyDensity(pt); }
    
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() {
        return new FEFiberMaterialPoint(m_pFint->CreateMaterialPointData());
    }
    
public:
    FEPropertyT<FEElasticFiberMaterialUC>   m_pFmat;    // pointer to fiber material
	FEPropertyT<FEFiberDensityDistribution> m_pFDD;     // pointer to fiber density distribution
	FEPropertyT<FEFiberIntegrationSchemeUC> m_pFint;    // pointer to fiber integration scheme
};

#endif /* defined(__FEBioMech__FEContinuousFiberDistributionUC__) */
