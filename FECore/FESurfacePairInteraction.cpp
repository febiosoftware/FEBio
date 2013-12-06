#include "stdafx.h"
#include "FESurfacePairInteraction.h"

//-----------------------------------------------------------------------------
FESurfacePairInteraction::FESurfacePairInteraction(FEModel* pfem) : FEModelComponent(FESURFACEPAIRINTERACTION_ID, pfem)
{
	m_nID = -1;
}
