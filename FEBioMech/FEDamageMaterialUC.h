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
#include "FEDamageCriterion.h"
#include "FEDamageCDF.h"

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
	bool Init();
    
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();
    
    // get the elastic material
    FEUncoupledMaterial* GetElasticMaterial() { return m_pBase; }
    
public:   
	//! Set the local coordinate system for a material point (overridden from FEMaterial)
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp);
    
public:
    FEPropertyT<FEUncoupledMaterial>    m_pBase;    // base elastic material
	FEPropertyT<FEDamageCDF>            m_pDamg;    // damage model
	FEPropertyT<FEDamageCriterion>      m_pCrit;    // damage criterion

	DECLARE_FECORE_CLASS();
};

#endif /* defined(__FEBioMech__FEDamageMaterialUC__) */
