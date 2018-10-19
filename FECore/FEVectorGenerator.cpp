#include "stdafx.h"
#include "FEVectorGenerator.h"
#include "FEMaterialPoint.h"
#include "FEElement.h"
#include "FEDomain.h"
#include "FENode.h"
#include "FEDomainMap.h"

FEVectorGenerator::FEVectorGenerator(FEModel* fem) : FECoreBase(fem, FEVECTORGENERATOR_ID)
{

}

FEVectorGenerator::~FEVectorGenerator()
{

}

//=================================================================================================
BEGIN_FECORE_CLASS(FELocalVectorGenerator, FEVectorGenerator)
	ADD_PARAMETER(m_n, 2, "local");
END_FECORE_CLASS();

FELocalVectorGenerator::FELocalVectorGenerator(FEModel* fem) : FEVectorGenerator(fem)
{
	m_n[0] = m_n[1] = 0;
}

bool FELocalVectorGenerator::Init()
{
	if ((m_n[0] == 0) && (m_n[1] == 0)) m_n[1] = 1;
	return FEVectorGenerator::Init();
}

vec3d FELocalVectorGenerator::GetVector(const FEMaterialPoint& mp)
{
	FEElement* el = mp.m_elem; assert(el);

	FEMeshPartition* dom = el->GetMeshPartition();
	vec3d r0 = dom->Node(el->m_lnode[m_n[0]]).m_r0;
	vec3d r1 = dom->Node(el->m_lnode[m_n[1]]).m_r0;

	vec3d n = r1 - r0;
	n.unit();

	return n;
}

//=================================================================================================
BEGIN_FECORE_CLASS(FEConstVectorGenerator, FEVectorGenerator)
	ADD_PARAMETER(m_vec, "vector");
END_FECORE_CLASS();

FEConstVectorGenerator::FEConstVectorGenerator(FEModel* fem) : FEVectorGenerator(fem)
{
	m_vec = vec3d(1, 0, 0);
}

bool FEConstVectorGenerator::Init()
{
	m_vec.unit();
	return true;
}

vec3d FEConstVectorGenerator::GetVector(const FEMaterialPoint& mp)
{
	return m_vec;
}

//=================================================================================================
BEGIN_FECORE_CLASS(FESphericalVectorGenerator, FEVectorGenerator)
	ADD_PARAMETER(m_center, "center");
	ADD_PARAMETER(m_vector, "vector");
END_FECORE_CLASS();

FESphericalVectorGenerator::FESphericalVectorGenerator(FEModel* fem) : FEVectorGenerator(fem)
{
	m_center = vec3d(0, 0, 0);
	m_vector = vec3d(1, 0, 0);
}

bool FESphericalVectorGenerator::Init()
{
	// Make sure the vector is a unit vector
	m_vector.unit();
	return true;
}

vec3d FESphericalVectorGenerator::GetVector(const FEMaterialPoint& mp)
{
	vec3d a = mp.m_r0 - m_center;
	a.unit();

	// setup the rotation
	vec3d e1(1, 0, 0);
	quatd q(e1, a);

	vec3d v = m_vector;
//	v.unit();	
	q.RotateVector(v);

	return v;
}

//=================================================================================================
BEGIN_FECORE_CLASS(FECylindricalVectorGenerator, FEVectorGenerator)
	ADD_PARAMETER(m_center, "center");
	ADD_PARAMETER(m_axis  , "axis"  )
	ADD_PARAMETER(m_vector, "vector");
END_FECORE_CLASS();

FECylindricalVectorGenerator::FECylindricalVectorGenerator(FEModel* fem) : FEVectorGenerator(fem)
{
	m_center = vec3d(0, 0, 0);
	m_axis   = vec3d(0, 0, 1);
	m_vector = vec3d(1, 0, 0);
}

bool FECylindricalVectorGenerator::Init()
{
	// Make sure the axis and vector are unit vectors
	m_axis.unit();
	m_vector.unit();
	return true;
}

vec3d FECylindricalVectorGenerator::GetVector(const FEMaterialPoint& mp)
{
	vec3d p = mp.m_r0 - m_center;

	// find the vector to the axis
	vec3d b = p - m_axis*(m_axis*p);
	b.unit();

	// setup the rotation
	vec3d e1(1, 0, 0);
	quatd q(e1, b);

	vec3d r = m_vector;
//	r.unit();	
	q.RotateVector(r);

	return r;
}

//=================================================================================================
FEUserVectorGenerator::FEUserVectorGenerator(FEModel* fem) : FEVectorGenerator(fem), m_data(nullptr)
{

}

bool FEUserVectorGenerator::Init()
{
	if (m_data == nullptr) return false;
	return FEVectorGenerator::Init();
}

void FEUserVectorGenerator::SetData(FEDomainMap* map)
{
	m_data = map;
}

vec3d FEUserVectorGenerator::GetVector(const FEMaterialPoint& mp)
{
	vec3d n = m_data->valueVec3d(mp);
	n.unit();
	return n;
}
