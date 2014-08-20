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
    FEContinuousFiberDistributionUC(FEModel* pfem) : FEUncoupledMaterial(pfem) {}
    ~FEContinuousFiberDistributionUC() {}
    
    // Initialization
    void Init();
    
	// get a material parameter
	FEParam* GetParameter(const ParamString& s);
    
public:
	//! get the number of material properties
	int Properties();
    
	//! get a specific material property
	FECoreBase* GetProperty(int i);
    
	//! find a material property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);
    
	//! set a material property (returns false on error)
	bool SetProperty(int i, FECoreBase* pm);
    
public:
	//! calculate stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) { return m_pFint->Stress(pt); }
    
	//! calculate tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) { return m_pFint->Tangent(pt); }
    
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() {
        return new FEFiberMaterialPoint(m_pFint->CreateMaterialPointData());
    }
    
public:
    FEElasticFiberMaterialUC*   m_pFmat;    // pointer to fiber material
    FEFiberDensityDistribution* m_pFDD;     // pointer to fiber density distribution
    FEFiberIntegrationSchemeUC* m_pFint;    // pointer to fiber integration scheme
};

#endif /* defined(__FEBioMech__FEContinuousFiberDistributionUC__) */
