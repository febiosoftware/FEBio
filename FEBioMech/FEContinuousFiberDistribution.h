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
    FEContinuousFiberDistribution(FEModel* pfem) : FEElasticMaterial(pfem) {
        m_a = vec3d(1,0,0); m_d = vec3d(0,1,0); m_Q = mat3dd(1);
    }
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
    
    // return local material axes
    mat3d LocalMatAxes() { return m_Q; }
    
private:
    vec3d   m_a;        // material axes relative to local axes
    vec3d   m_d;
    mat3d   m_Q;        // local orientation
    
public:
    FEElasticFiberMaterial*     m_pFmat;    // pointer to fiber material
    FEFiberDensityDistribution* m_pFDD;     // pointer to fiber density distribution
    FEFiberIntegrationScheme*   m_pFint;    // pointer to fiber integration scheme
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
