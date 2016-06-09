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
	void Init();

	//! Update material point data
	void Update(const FETimeInfo& timeInfo);

	//! Serialize data to archive
	void Serialize(DumpStream& ar);

public:
	mat3ds	m_se;	//!< elastic Cauchy stress
	mat3ds	m_Sep;	//!< elastic 2nd PK stress at previous time

	mat3ds	m_H[MAX_TERMS];		//!< internal variables
	mat3ds	m_Hp[MAX_TERMS];	//!< internal variables at previous timestep
    
//	double	m_sed;	//!< elastic strain energy density
//	double	m_sedp;	//!< elastic strain energy density at previous time
    
//	double	m_Hsed[MAX_TERMS];	//!< sed internal variables
//	double	m_Hsedp[MAX_TERMS];	//!< sed internal variables at previous timestep
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

	//! get the elastic base material
	FEElasticMaterial* GetBaseMaterial();

	//! Set the base material
	void SetBaseMaterial(FEElasticMaterial* pbase);

	//! Set the local coordinate system for a material point (overridden from FEMaterial)
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp);

public:
	//! stress function
	mat3ds Stress(FEMaterialPoint& pt);

	//! tangent function
	tens4ds Tangent(FEMaterialPoint& pt);

	//! strain energy density
	double StrainEnergyDensity(FEMaterialPoint& pt);
    
    // returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();

public: 
	// material parameters
	double	m_g0;			//!< intitial visco-elastic coefficient
	double	m_g[MAX_TERMS];	//!< visco-elastic coefficients
	double	m_t[MAX_TERMS];	//!< relaxation times

private:
	FEPropertyT<FEElasticMaterial>	m_Base;	//!< pointer to elastic solid material

public:
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
