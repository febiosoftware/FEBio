#include "stdafx.h"
#include "FESurfacePairInteraction.h"

int FESurfacePairInteraction::m_ncount = 0;

//-----------------------------------------------------------------------------
FESurfacePairInteraction::FESurfacePairInteraction(FEModel* pfem) : FEModelComponent(FESURFACEPAIRINTERACTION_ID, pfem)
{
	m_nID = m_ncount++;
}
