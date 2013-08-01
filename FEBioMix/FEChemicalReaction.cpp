#include "stdafx.h"
#include "FEChemicalReaction.h"

//-----------------------------------------------------------------------------
FEChemicalReaction::FEChemicalReaction()
{
	m_Vovr = false; 
	m_pMP = 0; 
}

//-----------------------------------------------------------------------------
void FEChemicalReaction::Init() 
{
	if (m_pFwd) InitializeReactionRate(m_pFwd);
	if (m_pRev) InitializeReactionRate(m_pRev);
}

//-----------------------------------------------------------------------------
// initialize chemical reaction rate
void FEChemicalReaction::InitializeReactionRate(FEReactionRate* m_pRate)
{
	m_pRate->m_pReact = this; 
}
