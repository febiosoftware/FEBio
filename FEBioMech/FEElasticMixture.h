#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material point data for mixtures
//!
class FEElasticMixtureMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEElasticMixtureMaterialPoint();

	//! Copy material point data
	FEMaterialPoint* Copy();

	//! material point initialization
	void Init(bool bflag);

	//! data serialization
	void Serialize(DumpFile& ar);

	//! stream material point data
	void ShallowCopy(DumpStream& dmp, bool bsave);

public:
	vector<double>				m_w;	//!< material weights
	vector<FEMaterialPoint*>	m_mp;	//!< material point data for mixture components
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
	FEMaterialPoint* CreateMaterialPointData();

	// return number of materials
	int Materials() const { return m_pMat.size(); }

	// return a material component
	FEElasticMaterial* GetMaterial(int i) { return m_pMat[i]; }

	// Add a material component
	void AddMaterial(FEElasticMaterial* pm);

	// get a material parameter
	FEParam* GetParameter(const ParamString& s);

public:
	//! return number of properties
	int Properties() { return (int) m_pMat.size(); }

	//! return a material property
	FEMaterial* GetProperty(int i) { return m_pMat[i]; }

	//! find a material property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! set a material property (returns false on error)
	bool SetProperty(int i, FEMaterial* pm);
	
public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);
		
	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);
		
	//! data initialization and checking
	void Init();

private:
	vector <FEElasticMaterial*>	m_pMat;	//!< pointers to elastic materials
};
