#include "stdafx.h"
#include "FETrussDomain.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
void FETrussDomain::Activate()
{
	for (int i=0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.m_bexclude == false)
		{
			if (node.m_rid < 0)
			{
				node.m_ID[DOF_X] = DOF_ACTIVE;
				node.m_ID[DOF_Y] = DOF_ACTIVE;
				node.m_ID[DOF_Z] = DOF_ACTIVE;
			}
		}
	}
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
