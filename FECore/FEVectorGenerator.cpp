#include "stdafx.h"
#include "FEVectorGenerator.h"
#include "FEMaterialPoint.h"
#include "FEElement.h"
#include "FEDomain.h"
#include "FENode.h"
#include "FEDomainMap.h"

FEVectorGenerator::FEVectorGenerator(FEModel* fem) : FECoreBase(FEVECTORGENERATOR_ID)
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

	FEDomain* dom = el->GetDomain();
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
