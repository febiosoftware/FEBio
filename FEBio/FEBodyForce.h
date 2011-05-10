#pragma once
#include "FECore/FEMaterialPoint.h"
#include "MathParser.h"
#include "FECore/DumpFile.h"
#include "FEMesh.h"

class FEM;

//-----------------------------------------------------------------------------
//! This class is the base class for body forces
//! Derived classes need to implement the force and stiffness functions.
//
class FEBodyForce
{
public:
	FEBodyForce()
	{
		s[0] = s[1] = s[2] = 0.0;
		lc[0] = -1; lc[1] = -1; lc[2] = -1;
	}

	//! calculate the body force at a material point
	virtual vec3d force(FEMaterialPoint& pt) = 0;
	virtual mat3ds stiffness(FEMaterialPoint& pt) = 0;

	virtual void Serialize(DumpFile& ar);

	virtual void Init(){}
	virtual void Update(){}

public:
	double	s[3];		// scale factor
	int		lc[3];		// load curve number
};

//-----------------------------------------------------------------------------
//! This class defines a deformation-independent constant force (e.g. gravity)

//! Note that the returned force is constanct. Use the scale factors and load
//! curves to define the intensity
class FEConstBodyForce : public FEBodyForce
{
public:
	vec3d force(FEMaterialPoint& pt) { return vec3d(1,1,1); }
	mat3ds stiffness(FEMaterialPoint& pt) { return mat3ds(0,0,0,0,0,0); }
};

//-----------------------------------------------------------------------------
//! This class defines a non-homogeneous force, i.e. the force depends
//! on the spatial position
class FENonConstBodyForce : public FEBodyForce
{
public:
	FENonConstBodyForce();
	vec3d force(FEMaterialPoint& pt);
	mat3ds stiffness(FEMaterialPoint& pt);
	void Serialize(DumpFile& ar);

public:
	char	m_sz[3][256];
};

//-----------------------------------------------------------------------------
//! This class defines a centrigufal force

class FECentrifugalBodyForce : public FEBodyForce
{
public:
	vec3d force(FEMaterialPoint& mp) {
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		mat3ds K = stiffness(mp);
		return K*(pt.rt - c);
	}
	mat3ds stiffness(FEMaterialPoint& mp) { return mat3dd(1) - dyad(n); }
	void Serialize(DumpFile& ar);
	
public:
	vec3d	n;	// rotation axis
	vec3d	c;	// point on axis of rotation (e.g., center of rotation)
};

//-----------------------------------------------------------------------------
class FEPointBodyForce : public FEBodyForce
{
public:
	// type of force center:
	// POINT = a global point, rigid or not
	// NODE = a node of the mesh
	enum { POINT, NODE };

public:
	FEPointBodyForce(FEM* pfem);

	vec3d force(FEMaterialPoint& mp);
	mat3ds stiffness(FEMaterialPoint& mp);

	void Serialize(DumpFile& ar);

	void Init();
	void Update();

public:
	FEM*	m_pfem;
	double	m_a, m_b;
	vec3d	m_rc;
	int		m_rlc[3];
	int		m_ntype;
	int		m_inode;

	bool	m_brigid;

	FESolidElement* m_pel;	// element in which point m_r0 lies
	double			m_rs[3];	// isoparametric coordinates
};
