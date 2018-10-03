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
	FEUncoupledViscoElasticMaterial(FEModel* pfem);

	// get the elastic base material
	FEElasticMaterial* GetBaseMaterial() { return m_pBase; }

	// set the elastic base material
	void SetBaseMaterial(FEUncoupledMaterial* pbase) { m_pBase = pbase; }

	//! Set the local coordinate system for a material point (overridden from FEMaterial)
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp) override;

public:
	//! data initialization and checking
	bool Init() override;
	
	//! deviatoric stress function
	mat3ds DevStress(FEMaterialPoint& pt) override;
	
	//! deviatoric tangent function
	tens4ds DevTangent(FEMaterialPoint& pt) override;
	
	//! deviatoric strain energy density function
	double DevStrainEnergyDensity(FEMaterialPoint& pt) override;
	
	//! returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() override;
	
public:
	double	m_t[MAX_TERMS];	//!< relaxation times
	double	m_g0;			//!< intitial visco-elastic coefficient
	double	m_g[MAX_TERMS];	//!< visco-elastic coefficients
	
private:
	FEPropertyT<FEUncoupledMaterial>	m_pBase;	//!< pointer to elastic solid material
	bool					m_binit;	//!< initialization flag
	
public:
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
