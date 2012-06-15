#pragma once
#include "FEUncoupledMaterial.h"
#include "FEViscoElasticMaterial.h"

//-----------------------------------------------------------------------------
//! This class implements a large deformation uncoupled visco-elastic material

class FEUncoupledViscoElasticMaterial :	public FEUncoupledMaterial
{
public:
	// NOTE: make sure that this parameter is the 
	//       same as the MAX_TERMS in the FEViscoElasticMaterialPoint class
	enum { MAX_TERMS = FEViscoElasticMaterialPoint::MAX_TERMS };
	
public:
	//! default constructor
	FEUncoupledViscoElasticMaterial();
	
	//! data initialization and checking
	void Init();
	
	//! deviatoric stress function
	mat3ds DevStress(FEMaterialPoint& pt);
	
	//! deviatoric tangent function
	tens4ds DevTangent(FEMaterialPoint& pt);
	
	// returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() 
	{ 
		return new FEViscoElasticMaterialPoint(m_pBase->CreateMaterialPointData());
	}
	
	// get the elastic material
	FEElasticMaterial* GetElasticMaterial() { return m_pBase; }
	
public:
	double	m_t[MAX_TERMS];	//!< relaxation times
	double	m_g0;			//!< intitial visco-elastic coefficient
	double	m_g[MAX_TERMS];	//!< visco-elastic coefficients
	
public:
	FEUncoupledMaterial*	m_pBase;	//!< pointer to elastic solid material
	
public:
	// declare as registered
	DECLARE_REGISTERED(FEUncoupledViscoElasticMaterial);
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
