#pragma once
#include "FEBodyForce.h"
#include "FECore/FEElement.h"

//-----------------------------------------------------------------------------
class FEPointBodyForce : public FEBodyForce
{
public:
	FEPointBodyForce(FEModel* pfem);

	vec3d force(FEMaterialPoint& mp);
	mat3ds stiffness(FEMaterialPoint& mp);

	void Serialize(DumpStream& ar);

	bool Init();
	void Update();

public:
	double	m_a, m_b;
	vec3d	m_rc;
	
	int		m_inode;

	bool	m_brigid;

	FESolidElement* m_pel;		//!< element in which point m_r0 lies
	double			m_rs[3];	//!< isoparametric coordinates

	DECLARE_PARAMETER_LIST();
};
