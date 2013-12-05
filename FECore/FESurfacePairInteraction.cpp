#include "stdafx.h"
#include "FESurfacePairInteraction.h"

//-----------------------------------------------------------------------------
FESurfacePairInteraction::FESurfacePairInteraction(FEModel* pfem) : FECoreBase(FESURFACEPAIRINTERACTION_ID), m_pfem(pfem)
{
	m_nID = -1;
	m_bactive = true;
}

//-----------------------------------------------------------------------------
bool FESurfacePairInteraction::IsActive()
{ 
	return m_bactive; 
}

//-----------------------------------------------------------------------------
void FESurfacePairInteraction::Activate()
{ 
	m_bactive = true; 
}

//-----------------------------------------------------------------------------
void FESurfacePairInteraction::Deactivate()
{ 
	m_bactive = false; 
}
