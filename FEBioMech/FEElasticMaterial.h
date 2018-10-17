#pragma once
#include "FESolidMaterial.h"
#include <FECore/FEMaterialPointProperty.h>

//-----------------------------------------------------------------------------
//! This class defines material point data for elastic materials.
class FEElasticMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEElasticMaterialPoint();

	//! Initialize material point data
	void Init();

	//! create a shallow copy
	FEMaterialPoint* Copy();

	//! serialize material point data
	void Serialize(DumpStream& ar);

public:
	mat3ds Strain() const;
	mat3ds SmallStrain() const;

	mat3ds RightCauchyGreen() const;
	mat3ds LeftCauchyGreen () const;

	mat3ds DevRightCauchyGreen() const;
	mat3ds DevLeftCauchyGreen () const;
    
    mat3ds RateOfDeformation() const { return m_L.sym(); }

	mat3ds pull_back(const mat3ds& A) const;
	mat3ds push_forward(const mat3ds& A) const;

	tens4ds pull_back(const tens4ds& C) const;
	tens4ds push_forward(const tens4ds& C) const;

public:
    bool    m_buncoupled;   //!< set to true if this material point was created by an uncoupled material
    mat3d   m_Q;            //!< local material orientation
    
	// deformation data at intermediate time
    vec3d   m_rt;   //!< spatial position
	mat3d	m_F;	//!< deformation gradient
	double	m_J;	//!< determinant of F
    vec3d   m_v;    //!< velocity
    vec3d   m_a;    //!< acceleration
    mat3d   m_L;    //!< spatial velocity gradient

	// solid material data
	mat3ds		m_s;		//!< Cauchy stress
	mat3ds		m_s0;		//!< Initial stress (only used by linear solid solver)
    
    // current time data
    double      m_Wt;       //!< strain energy density at current time
    
    // previous time data
    double      m_Wp;       //!< strain energy density
};

//-----------------------------------------------------------------------------
class FEMatAxis : public FEMaterialPointProperty_T<FEElasticMaterialPoint, mat3d>
{
public:
	mat3d& data(FEElasticMaterialPoint& mp) override
	{
		return mp.m_Q;
	}
};

//-----------------------------------------------------------------------------
//! Base class for (hyper-)elastic materials

class FEElasticMaterial : public FESolidMaterial
{
public:
	//! constructor 
	FEElasticMaterial(FEModel* pfem);

	//! destructor
	~FEElasticMaterial();

	//! create material point data for this material
	virtual FEMaterialPoint* CreateMaterialPointData() { return new FEElasticMaterialPoint; }

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt);
    
	//! Get the elastic component
	FEElasticMaterial* GetElasticMaterial() { return this; }

	//! Set the local coordinate system for a material point (overridden from FEMaterial)
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp);

protected:
	FEMatAxis	m_Q;

	DECLARE_FECORE_CLASS();
};
