#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material point class for the PreStrainElastic material.
//! This class allows the user to define a pre-strain deformation gradient for
//! each integration point.
class FEPreStrainMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEPreStrainMaterialPoint(FEMaterialPoint* pt);

	//! Initialize
	void Init(bool bflag);

	//! create a shallow copy
	FEMaterialPoint* Copy();

	//! serialize material point data
	void Serialize(DumpFile& ar);

	//! stream material point data
	void ShallowCopy(DumpStream& dmp, bool bsave);

public:
	mat3d	m_Fp;	//!< prestrain deformation gradient
};

//-----------------------------------------------------------------------------
//! This material applies a user-defined prestrain deformation gradient
//! before evaluating the stress and tangents. 
class FEPreStrainElastic : public FEElasticMaterial
{
public:
	//! constructor
	FEPreStrainElastic(FEModel* pfem);

	//! Initialization
	void Init();

	//! Set the base material
	void SetBaseMaterial(FEElasticMaterial* pbase);

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();

public:
	//! return number of properties
	int Properties();

	//! return a material property
	FECoreBase* GetProperty(int i);

	//! find a material property index ( returns <0 for error)
	virtual int FindPropertyIndex(const char* szname);

	//! set a material property (returns false on error)
	virtual bool SetProperty(int i, FECoreBase* pm);

public:
	//! Cauchy stress 
	mat3ds Stress(FEMaterialPoint& mp);

	//! spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp);

private:
	FEElasticMaterial*	m_pmat;	//!< elastic base material

	mat3d	m_Fp; //!< pre-strain deformationg gradient

	DECLARE_PARAMETER_LIST();
};
