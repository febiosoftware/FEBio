#pragma once
#include "FEUncoupledMaterial.h"
#include "FEElasticMixture.h"

//-----------------------------------------------------------------------------
//! Uncoupled elastic mixtures

//! This class describes a mixture of uncoupled elastic solids.  The user must declare
//! uncoupled elastic solids that can be combined within this class.  The stress and
//! tangent tensors evaluated in this class represent the sum of the respective
//! tensors of all the solids forming the mixture.
//! \todo This class defines two accessor interfaces. Modify to use the FEMaterial interface only.

class FEUncoupledElasticMixture : public FEUncoupledMaterial
{
public:
	FEUncoupledElasticMixture(FEModel* pfem);

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();

	// return number of materials
	int Materials() { return (int)m_pMat.size(); }

	// return a material component
	FEUncoupledMaterial* GetMaterial(int i) { return m_pMat[i]; }

	// Add a material component
	void AddMaterial(FEUncoupledMaterial* pm);

	//! Set the local coordinate system for a material point (overridden from FEMaterial)
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp);

public:
	//! calculate stress at material point
	mat3ds DevStress(FEMaterialPoint& pt);
	
	//! calculate tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt);
	
	//! calculate strain energy density at material point
	double DevStrainEnergyDensity(FEMaterialPoint& pt);
    
	//! data initialization and checking
	bool Init();

private:
	FEVecPropertyT<FEUncoupledMaterial>	m_pMat;	//!< pointers to elastic materials

	DECLARE_FECORE_CLASS();
};
