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
int FEChemicalReaction::Properties()
{
	int n = 0;
	if (m_pFwd) n++;
	if (m_pRev) n++;
	return n;
}

//-----------------------------------------------------------------------------
FEMaterial* FEChemicalReaction::GetProperty(int i)
{
	int n = Properties();
	if (n == 1)
	{
		if (m_pFwd) return m_pFwd;
		if (m_pRev) return m_pRev;
	}
	else if (n==2)
	{
		if (i==0) return m_pFwd;
		if (i==1) return m_pRev;
	}
	return 0;
}

//-----------------------------------------------------------------------------
// initialize chemical reaction rate
void FEChemicalReaction::InitializeReactionRate(FEReactionRate* m_pRate)
{
	m_pRate->m_pReact = this; 
}
