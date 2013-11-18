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
	FEUncoupledElasticMixture(FEModel* pfem) : FEUncoupledMaterial(pfem) {}

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();

	// return number of materials
	int Materials() const { return m_pMat.size(); }

	// return a material component
	FEUncoupledMaterial* GetMaterial(int i) { return m_pMat[i]; }

	// Add a material component
	void AddMaterial(FEUncoupledMaterial* pm);

	// get a material parameter
	FEParam* GetParameter(const ParamString& s);

public:
	//! get number of material properties
	int Properties() { return m_pMat.size(); }

	//! return a material property
	FEMaterial* GetProperty(int n) { return m_pMat[n]; }

	//! find a material property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! set a material property (returns false on error)
	bool SetProperty(int i, FEMaterial* pm);
	
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
