#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material point data for mixtures
//!
class FEElasticMixtureMaterialPoint : public FEMaterialPointArray
{
public:
	//! constructor
	FEElasticMixtureMaterialPoint();

	//! Copy material point data
	FEMaterialPoint* Copy();

	//! material point initialization
	void Init();

	//! data serialization
	void Serialize(DumpStream& ar);

public:
	vector<double>				m_w;	//!< material weights
};

//-----------------------------------------------------------------------------
//! Elastic mixtures

//! This class describes a mixture of elastic solids.  The user must declare
//! elastic solids that can be combined within this class.  The stress and
//! tangent tensors evaluated in this class represent the sum of the respective
//! tensors of all the solids forming the mixture.

//! \todo This class defines two accessor interfaces. Modify to use the FEMaterial interface only.
class FEElasticMixture : public FEElasticMaterial
{
public:
	FEElasticMixture(FEModel* pfem);

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() override;

	// return number of materials
	int Materials() { return (int)m_pMat.size(); }

	// return a material component
	FEElasticMaterial* GetMaterial(int i) { return m_pMat[i]; }

	// Add a material component
	void AddMaterial(FEElasticMaterial* pm);

    //! specialized material points
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override;
    
	//! Set the local coordinate system for a material point (overridden from FEMaterial)
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp) override;

public:
   
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) override;
		
	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;
		
	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
private:
	FEVecPropertyT<FEElasticMaterial>	m_pMat;	//!< pointers to elastic materials
};
