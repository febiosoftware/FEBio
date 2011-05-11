#include "FECore/FEBodyForce.h"
#include "FECore/FEBioFactory.h"
#include "FECore/febio.h"

//-----------------------------------------------------------------------------
class FEPointBodyForce : public FEBodyForce
{
public:
	// type of force center:
	// POINT = a global point, rigid or not
	// NODE = a node of the mesh
	enum { POINT, NODE };

public:
	FEPointBodyForce(FEModel* pfem);

	vec3d force(FEMaterialPoint& mp);
	mat3ds stiffness(FEMaterialPoint& mp);

	void Serialize(DumpFile& ar);

	void Init();
	void Update();

public:
	double	m_a, m_b;
	vec3d	m_rc;
	int		m_rlc[3];
	int		m_ntype;
	int		m_inode;

	bool	m_brigid;

	FESolidElement* m_pel;		//!< element in which point m_r0 lies
	double			m_rs[3];	//!< isoparametric coordinates

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
class FEPointBodyForceFactory : public FEBodyForceFactory
{
public:
	FEBodyForce* Create(FEModel* pfem) { return new FEPointBodyForce(pfem); }
};
