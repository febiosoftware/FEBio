#pragma once
#include "FECore/FEBodyForce.h"
#include "MathParser.h"

//-----------------------------------------------------------------------------
//! This class defines a deformation-independent constant force (e.g. gravity)

//! Note that the returned force is constanct. Use the scale factors and load
//! curves to define the intensity
class FEConstBodyForce : public FEBodyForce
{
public:
	FEConstBodyForce(FEModel* pfem) : FEBodyForce(pfem) {}
	vec3d force(FEMaterialPoint& pt) { return vec3d(1,1,1); }
	mat3ds stiffness(FEMaterialPoint& pt) { return mat3ds(0,0,0,0,0,0); }
};

//-----------------------------------------------------------------------------
//! This class defines a non-homogeneous force, i.e. the force depends
//! on the spatial position
class FENonConstBodyForce : public FEBodyForce
{
public:
	FENonConstBodyForce(FEModel* pfem);
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
	FECentrifugalBodyForce(FEModel* pfem) : FEBodyForce(pfem){}
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
