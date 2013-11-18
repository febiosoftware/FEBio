#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material point data for visco-elastic materials
class FEViscoElasticMaterialPoint : public FEMaterialPoint
{
public:
	enum { MAX_TERMS = 6 };

public:
	//! constructor
	FEViscoElasticMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) {}

	//! copy material point data
	FEMaterialPoint* Copy();

	//! Initialize material point data
	void Init(bool bflag);

	//! Serialize data to archive
	void Serialize(DumpFile& ar);

public:
	mat3ds	m_se;	//!< elastic Cauchy stress
	mat3ds	m_Sep;	//!< elastic 2nd PK stress at previous time

	mat3ds	m_H[MAX_TERMS];		//!< internal variables
	mat3ds	m_Hp[MAX_TERMS];	//!< internal variables at previous timestep
};


//-----------------------------------------------------------------------------
//! This class implements a large deformation visco-elastic material
//
class FEViscoElasticMaterial :	public FEElasticMaterial
{
public:
	// NOTE: make sure that this parameter is the 
	//       same as the MAX_TERMS in the FEViscoElasticMaterialPoint class
	enum { MAX_TERMS = FEViscoElasticMaterialPoint::MAX_TERMS };

public:
	//! default constructor
	FEViscoElasticMaterial(FEModel* pfem);

	//! Get a parameter
	FEParam* GetParameter(const ParamString& s);

	//! get the elastic base material \todo I want to call this GetElasticMaterial, but this name is being used
	FEElasticMaterial* GetBaseMaterial() { return m_pBase; }

	//! Set the base material
	void SetBaseMaterial(FEElasticMaterial* pbase) { m_pBase = pbase; }

public:
	//! return number of properties
	int Properties();

	//! return a material property
	FEMaterial* GetProperty(int i);

	//! find a material property index ( returns <0 for error)
	virtual int FindPropertyIndex(const char* szname);

	//! set a material property (returns false on error)
	virtual bool SetProperty(int i, FEMaterial* pm);

public:
	//! data initialization
	void Init();

	//! stress function
	mat3ds Stress(FEMaterialPoint& pt);

	//! tangent function
	tens4ds Tangent(FEMaterialPoint& pt);

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();

public: 
	// material parameters
	double	m_g0;			//!< intitial visco-elastic coefficient
	double	m_g[MAX_TERMS];	//!< visco-elastic coefficients
	double	m_t[MAX_TERMS];	//!< relaxation times

private:
	FEElasticMaterial*	m_pBase;	//!< pointer to elastic solid material

public:
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
