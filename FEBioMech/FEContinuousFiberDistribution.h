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
    FEContinuousFiberDistribution(FEModel* pfem) : FEElasticMaterial(pfem) {}
    ~FEContinuousFiberDistribution() {}
    
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
	mat3ds Stress(FEMaterialPoint& pt) { return m_pFint->Stress(pt); }
    
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) { return m_pFint->Tangent(pt); }
    
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() {
        return new FEFiberMaterialPoint(m_pFint->CreateMaterialPointData());
    }
    
public:
    FEElasticFiberMaterial*     m_pFmat;    // pointer to fiber material
    FEFiberDensityDistribution* m_pFDD;     // pointer to fiber density distribution
    FEFiberIntegrationScheme*   m_pFint;    // pointer to fiber integration scheme
};
