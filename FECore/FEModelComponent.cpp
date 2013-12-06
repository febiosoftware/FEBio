#include "stdafx.h"
#include "FEModelComponent.h"

//-----------------------------------------------------------------------------
FEModelComponent::FEModelComponent(SUPER_CLASS_ID sid, FEModel* pfem) : FECoreBase(sid)
{
	m_pfem = pfem;
	m_bactive = true;
}

//-----------------------------------------------------------------------------
FEModelComponent::~FEModelComponent()
{

}

//-----------------------------------------------------------------------------
FEModel* FEModelComponent::GetFEModel()
{
	return m_pfem;
}

//-----------------------------------------------------------------------------
bool FEModelComponent::IsActive()
{ 
	return m_bactive; 
}

//-----------------------------------------------------------------------------
void FEModelComponent::Activate()
{ 
	m_bactive = true; 
}

//-----------------------------------------------------------------------------
void FEModelComponent::Deactivate()
{ 
	m_bactive = false; 
}
