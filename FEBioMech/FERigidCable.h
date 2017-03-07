#pragma once
#include <FECore/FEModelLoad.h>

class FERigidBody;

class FERigidCable : public FEModelLoad
{
	class FECablePoint : public FECoreBase
	{
	public:
		FECablePoint(FEModel* fem) : FECoreBase(FEOBJECT_ID){}

		bool SetAttribute(const char* szatt, const char* szval);

	public:
		int		m_rb;	//!< rigid body ID
		vec3d	m_pos;	//!< position of attachment point

		DECLARE_PARAMETER_LIST();
	};

public:
	FERigidCable(FEModel* fem);

	//! initialization
	bool Init();

	//! Residual
	void Residual(FEGlobalVector& R, const FETimeInfo& tp);

	//! Stiffness matrix
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp);

	//! override for building points list
	FECoreBase* GetProperty(int n);

private:
	void applyRigidForce(FERigidBody& rb, const vec3d& F, const vec3d& d, FEGlobalVector& R);

private:
	double	m_force;		//!< magnitude of force (i.e. tension in cable)
	vec3d	m_forceDir;		//!< direction of force at cable's end
	bool	m_brelative;	//!< positions are defined relative w.r.t. rigid body's COM or not
	FEVecPropertyT<FECablePoint>	m_points;

private:
	DECLARE_PARAMETER_LIST();
};
