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
	FEUncoupledElasticMixture() {}

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();

	// return number of materials
	int Materials() const { return m_pMat.size(); }

	// return a material component
	FEUncoupledMaterial* GetMaterial(int i) { return m_pMat[i]; }

	// Add a material component
	void AddMaterial(FEUncoupledMaterial* pm) { m_pMat.push_back(pm); }

	// get a material parameter
	FEParam* GetParameter(const ParamString& s);

public:
	//! get number of material properties
	int Properties() { return m_pMat.size(); }

	//! return a material property
	FEMaterial* GetProperty(int n) { return m_pMat[n]; }
	
public:
	//! calculate stress at material point
	mat3ds DevStress(FEMaterialPoint& pt);
	
	//! calculate tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt);
	
	//! data initialization and checking
	void Init();

private:
	vector <FEUncoupledMaterial*>	m_pMat;	//!< pointers to elastic materials
};
