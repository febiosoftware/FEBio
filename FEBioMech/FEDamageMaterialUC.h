#pragma once
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
	mat3ds DevStress(FEMaterialPoint& pt) override;
    
	//! calculate tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;
    
	//! calculate strain energy density at material point
	double DevStrainEnergyDensity(FEMaterialPoint& pt) override;
    
    //! damage
    double Damage(FEMaterialPoint& pt);
    
	//! data initialization and checking
	bool Init() override;
    
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() override;
    
    // get the elastic material
    FEUncoupledMaterial* GetElasticMaterial() { return m_pBase; }
    
public:
    FEUncoupledMaterial*    m_pBase;    // base elastic material
	FEDamageCDF*            m_pDamg;    // damage model
	FEDamageCriterion*      m_pCrit;    // damage criterion

	DECLARE_FECORE_CLASS();
};
