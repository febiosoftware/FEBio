//
//  FEDamageMaterialUC.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/19/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEDamageMaterialUC__
#define __FEBioMech__FEDamageMaterialUC__

#include "FEUncoupledMaterial.h"
#include "FEDamageMaterialPoint.h"

class FEDamageCriterionUC;
class FEDamageCDF;

//-----------------------------------------------------------------------------
// This material models damage in any hyper-elastic materials.

class FEDamageMaterialUC : public FEUncoupledMaterial
{
public:
	FEDamageMaterialUC(FEModel* pfem);
    
public:
	//! calculate stress at material point
	mat3ds DevStress(FEMaterialPoint& pt);
    
	//! calculate tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt);
    
	//! calculate strain energy density at material point
	double DevStrainEnergyDensity(FEMaterialPoint& pt);
    
    //! damage
    double Damage(FEMaterialPoint& pt);
    
	//! data initialization and checking
	void Init();
    
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData()
	{
        return new FEDamageMaterialPoint(new FEElasticMaterialPoint, m_pCrit);
	}
    
    // get the elastic material
    FEUncoupledMaterial* GetElasticMaterial() { return m_pBase; }
    
public:
	// get a material parameter
	FEParam* GetParameter(const ParamString& s);
    
	//! get the number of material properties
	int Properties();
    
	//! get a specific material property
	FECoreBase* GetProperty(int i);
    
	//! find a material property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);
    
	//! set a material property (returns false on error)
	bool SetProperty(int i, FECoreBase* pm);
    
	//! Set the local coordinate system for a material point (overridden from FEMaterial)
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp);
    
	//! data serialization
	void Serialize(DumpFile& ar);

public:
    FEUncoupledMaterial*    m_pBase;    // base elastic material
    FEDamageCDF*            m_pDamg;    // damage model
    FEDamageCriterionUC*    m_pCrit;    // damage criterion
};

#endif /* defined(__FEBioMech__FEDamageMaterialUC__) */
