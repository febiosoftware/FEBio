#pragma once
#include "FEMaterial.h"
#include "FECoordSysMap.h"

//-----------------------------------------------------------------------------
//! This class defines material point data for elastic materials.
class FEElasticMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEElasticMaterialPoint();

	//! Initialize material point data
	void Init(bool bflag);

	//! create a shallow copy
	FEMaterialPoint* Copy();

	//! serialize material point data
	void Serialize(DumpFile& ar);

public:
	mat3ds Strain();
	mat3ds SmallStrain();

	mat3ds RightCauchyGreen();
	mat3ds LeftCauchyGreen ();

	mat3ds DevRightCauchyGreen();
	mat3ds DevLeftCauchyGreen ();

	mat3ds pull_back(const mat3ds& A);
	mat3ds push_forward(const mat3ds& A);

public:
	// position 
	vec3d	m_r0;	//!< material position
	vec3d	m_rt;	//!< spatial position

	// deformation data
	mat3d	m_F;	//!< deformation gradient
	double	m_J;	//!< determinant of F
	mat3d	m_Q;	//!< local material orientation

	// solid material data
	mat3ds		m_s;		//!< Cauchy stress
	mat3ds		m_s0;		//!< Initial stress (only used by linear solid solver)
	double		m_sed;		//!< strain energy density	\todo Is this a good place for this?
	double		m_rhor;		//!< current referential mass density
};

//-----------------------------------------------------------------------------
//! Base class for (hyper-)elastic materials

class FEElasticMaterial : public FESolidMaterial
{
public:
	//! constructor 
	FEElasticMaterial();

	//! destructor
	~FEElasticMaterial();

	//! Initialization
	void Init();

	//! Serialization
	void Serialize(DumpFile& ar);

	//! create material point data for this material
	virtual FEMaterialPoint* CreateMaterialPointData() { return new FEElasticMaterialPoint; }

	//! Get the elastic component
	FEElasticMaterial* GetElasticMaterial() { return this; }

public:
	bool			m_unstable;		//!< flag indicating whether material is unstable on its own
	FEMaterial*		m_pParent;		//!< pointer to parent	\todo This has to go!
	FECoordSysMap*	m_pmap;			//!< local material coordinate system

	DECLARE_PARAMETER_LIST();
};
