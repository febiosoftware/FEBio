#pragma once
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! This class implements a large deformation visco-elastic material

class FEViscoElasticMaterial :	public FEMaterial
{
public:
	// NOTE: make sure that this parameter is the 
	//       same as the MAX_TERMS in the FEViscoElasticMaterialPoint class
	enum { MAX_TERMS = FEViscoElasticMaterialPoint::MAX_TERMS };

public:
	//! default constructor
	FEViscoElasticMaterial();

	//! stress function
	mat3ds Stress(FEMaterialPoint& pt);

	//! tangent function
	void Tangent(double D[6][6], FEMaterialPoint& pt);

	// returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() 
	{ 
		return new FEViscoElasticMaterialPoint(m_pemat->CreateMaterialPointData());
	}

	//! return elastic component's bulk modulus
	double BulkModulus() { return m_pemat->BulkModulus(); }

	//! return elastic component's density
	double Density() { return m_pemat->Density(); }


public:
	int	m_nSolidMat;				//<! material id of elastic material
	FEElasticMaterial*	m_pemat;	//!< pointer to elastic material

	double	m_t[MAX_TERMS];	//!< relaxation times
	double	m_g0;			//!< intitial visco-elastic coefficient
	double	m_g[MAX_TERMS];	//!< visco-elastic coefficients

public:
	// declare as registered
	DECLARE_REGISTERED(FEViscoElasticMaterial);

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
