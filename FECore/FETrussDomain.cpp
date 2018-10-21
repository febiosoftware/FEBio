#include "stdafx.h"
#include "FETrussDomain.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
FETrussDomain::FETrussDomain(FEModel* fem) : FEDomain(FE_DOMAIN_TRUSS, fem)
{
}

//-----------------------------------------------------------------------------
void FETrussDomain::Create(int nsize, int elemType)
{
	m_Elem.resize(nsize);
	for (int i = 0; i < nsize; ++i)
	{
		FETrussElement& el = m_Elem[i];
		el.SetLocalID(i);
		el.SetMeshPartition(this);
	}

	if (elemType != -1)
		for (int i=0; i<nsize; ++i) m_Elem[i].SetType(elemType);
}

//-----------------------------------------------------------------------------
vec3d FETrussDomain::TrussNormal(FETrussElement& el)
{
	vec3d rt[2];
	rt[0] = m_pMesh->Node(el.m_node[0]).m_rt;
	rt[1] = m_pMesh->Node(el.m_node[1]).m_rt;

	vec3d a = rt[0];
	vec3d b = rt[1];
	vec3d n = b - a;
	n.unit();
	return n;
}

//-----------------------------------------------------------------------------
void FETrussDomain::ForEachTrussElement(std::function<void(FETrussElement& el)> f)
{
	int N = Elements();
	for (int i = 0; i < N; ++i)
	{
		FETrussElement& el = m_Elem[i];
		f(el);
	}
}
