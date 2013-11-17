#pragma once
#include "FEUncoupledMaterial.h"
#include "FEViscoElasticMaterial.h"

//-----------------------------------------------------------------------------
//! This class implements a large deformation uncoupled visco-elastic material
//
class FEUncoupledViscoElasticMaterial :	public FEUncoupledMaterial
{
public:
	// NOTE: make sure that this parameter is the 
	//       same as the MAX_TERMS in the FEViscoElasticMaterialPoint class
	enum { MAX_TERMS = FEViscoElasticMaterialPoint::MAX_TERMS };
	
public:
	//! default constructor
	FEUncoupledViscoElasticMaterial();

	// get the elastic material \todo Is this right to overload this function?
	FEElasticMaterial* GetElasticMaterial() { return m_pBase; }

	// get the elastic base material
	FEElasticMaterial* GetBaseMaterial() { return m_pBase; }

	// set the elastic base material
	void SetBaseMaterial(FEUncoupledMaterial* pbase) { m_pBase = pbase; }

	//! Find a material parameter
	FEParam* GetParameter(const ParamString& s);

public:
	//! return number of properties
	int Properties();

	//! return a material property
	FEMaterial* GetProperty(int i);

	//! find a material property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! set a material property (returns false on error)
	bool SetProperty(int i, FEMaterial* pm);

public:
	//! data initialization and checking
	void Init();
	
	//! deviatoric stress function
	mat3ds DevStress(FEMaterialPoint& pt);
	
	//! deviatoric tangent function
	tens4ds DevTangent(FEMaterialPoint& pt);
	
	//! returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();
	
public:
	double	m_t[MAX_TERMS];	//!< relaxation times
	double	m_g0;			//!< intitial visco-elastic coefficient
	double	m_g[MAX_TERMS];	//!< visco-elastic coefficients
	
private:
	FEUncoupledMaterial*	m_pBase;	//!< pointer to elastic solid material
	bool					m_binit;	//!< initialization flag
	
public:
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
