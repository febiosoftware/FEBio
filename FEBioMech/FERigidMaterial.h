#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Rigd body material data

//! Since rigid elements are skipped during the stiffness and residual calculations
//! we don't implement the Stress and Tangent functions
//! \todo make the m_rc a parameter
//! \todo Can I remove the m_bc variable?

class FERigidMaterial : public FESolidMaterial
{
public:
	FERigidMaterial(FEModel* pfem);

public:
	double	m_E;		//!< Young's modulus
	double	m_v;		//!< Poisson's ratio
	int		m_pmid;		//!< parent material ID

public:
	int		m_com;	//!< center of mass input flag
	vec3d	m_rc;	//!< center of mass
	int		m_nRB;	//!< rigid body ID (TODO: rigid materials can be assigned to mulitple rigid bodies, so does it make sense to store this?)

public:
	// inherited from FEMaterial
	bool IsRigid() const override { return true; }

	// override this function to set the COM logic
	void SetParameter(FEParam& p) override;

public:
	//! get the ID of the rigid body this material is assigned to (-1 if not)
	int GetRigidBodyID() { return m_nRB; }

	//! Set the rigid body ID this material is assigned to
	void SetRigidBodyID(int rid) { m_nRB = rid; }

public:
	//! Create a rigid material point
	FEMaterialPoint* CreateMaterialPointData() override { return new FEElasticMaterialPoint(); }

	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) override { return mat3ds(); }

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override { return tens4ds(); }

	//! data initialization
	bool Init() override;

	//! serialization
	void Serialize(DumpStream& ar) override;

	// declare a parameter list
	DECLARE_FECORE_CLASS();

private:
	bool	m_binit;	//!< flag for first initialization
};
