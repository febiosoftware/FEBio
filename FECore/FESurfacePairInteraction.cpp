#include "stdafx.h"
#include "FESurfacePairInteraction.h"

//-----------------------------------------------------------------------------
FESurfacePairInteraction::FESurfacePairInteraction(FEModel* pfem) : m_pfem(pfem)
{
	m_ntype = 0;
	m_nID = -1;
	m_bactive = true;
}

//-----------------------------------------------------------------------------
int FESurfacePairInteraction::Type()
{ 
	return m_ntype; 
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
