#include "stdafx.h"
#include "FEModelComponent.h"

//-----------------------------------------------------------------------------
//! The constructor takes two arguments: the SUPER_CLASS_ID which defines the 
//! type of the model component and a pointer to the FEModel object this component
//! belongs to.
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
