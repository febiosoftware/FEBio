#pragma once
#include "FECore/FEBodyForce.h"
#include "FECore/MathParser.h"
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! This class defines a deformation-independent constant force (e.g. gravity)
class FEConstBodyForce : public FEBodyForce
{
public:
	FEConstBodyForce(FEModel* pfem) : FEBodyForce(pfem) { m_f = vec3d(0,0,0); }
	vec3d force(FEMaterialPoint& pt) { return m_f; }
	mat3ds stiffness(FEMaterialPoint& pt) { return mat3ds(0,0,0,0,0,0); }

protected:
	vec3d	m_f;

	DECLARE_PARAMETER_LIST();
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

	DECLARE_PARAMETER_LIST();
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
		return K*((pt.m_rt - c)*(w*w));
	}
	mat3ds stiffness(FEMaterialPoint& mp) { return (mat3dd(1) - dyad(n))*(-w*w); }
	void Serialize(DumpFile& ar);
	
public:
	vec3d	n;	// rotation axis
	vec3d	c;	// point on axis of rotation (e.g., center of rotation)
	double	w;	// angular speed

	DECLARE_PARAMETER_LIST();
};
