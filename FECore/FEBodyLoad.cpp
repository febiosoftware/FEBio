#include "stdafx.h"
#include "FEBodyLoad.h"
#include "FEModel.h"
#include "FEModelParam.h"

//-----------------------------------------------------------------------------
FEBodyLoad::FEBodyLoad(FEModel* pfem) : FEModelComponent(FEBODYLOAD_ID, pfem)
{
}

//-----------------------------------------------------------------------------
FEBodyLoad::~FEBodyLoad()
{
}

//-----------------------------------------------------------------------------
//! initialization
bool FEBodyLoad::Init()
{
	// If the domain list is empty, add all the domains
	if (m_dom.empty())
	{
		FEMesh& mesh = GetFEModel()->GetMesh();
		for (int i=0; i<mesh.Domains(); ++i)
		{
			FEDomain* dom = &mesh.Domain(i);
			m_dom.push_back(dom);
		}
	}
	return FEModelComponent::Init(); 
}

//-----------------------------------------------------------------------------
//! update
void FEBodyLoad::Update()
{
}

//-----------------------------------------------------------------------------
int FEBodyLoad::Domains() const
{
	return (int) m_dom.size();
}

//-----------------------------------------------------------------------------
FEDomain* FEBodyLoad::Domain(int i)
{
	return m_dom[i];
}

//-----------------------------------------------------------------------------
void FEBodyLoad::AddDomain(FEDomain* dom)
{
	m_dom.push_back(dom);

	// add it do all the mapped parameters
	FEParameterList& PL = GetParameterList();
	FEParamIterator it = PL.first();
	for (int i = 0; i < PL.Parameters(); ++i, ++it)
	{
		FEParam& pi = *it;
		if (pi.type() == FE_PARAM_DOUBLE_MAPPED)
		{
			FEParamDouble& param = pi.value<FEParamDouble>();
			param.addDomain(dom);
		}
		if (pi.type() == FE_PARAM_VEC3D_MAPPED)
		{
			FEParamVec3& param = pi.value<FEParamVec3>();
			param.addDomain(dom);
		}
	}
}
