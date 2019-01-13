#pragma once
#include <FECore/FEMaterialPoint.h>
#include <FECore/tens4d.h>
#include <FECore/FEMaterialPointProperty.h>

//-----------------------------------------------------------------------------
//! This class defines material point data for elastic materials.
class FECORE_API FEElasticMaterialPoint : public FEMaterialPoint
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
    
	// deformation data at intermediate time
    vec3d   m_rt;   //!< spatial position
	mat3d	m_F;	//!< deformation gradient
	double	m_J;	//!< determinant of F
    vec3d   m_v;    //!< velocity
    vec3d   m_a;    //!< acceleration
    mat3d   m_L;    //!< spatial velocity gradient

	// solid material data
	mat3ds		m_s;		//!< Cauchy stress
    
    // current time data
    double      m_Wt;       //!< strain energy density at current time
    
    // previous time data
    double      m_Wp;       //!< strain energy density
};
