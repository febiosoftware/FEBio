#include "stdafx.h"
#include "FEChemicalReaction.h"
#include "FECore/FEElementTraits.h"
#include <stdlib.h>


//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEChemicalReaction, FEMaterial)
	ADD_PARAMETER(m_Vbar , FE_PARAM_DOUBLE, "Vbar");
	ADD_PARAMETER(m_vRtmp, FE_PARAM_INT   , "vR"  );
	ADD_PARAMETER(m_vPtmp, FE_PARAM_INT   , "vP"  );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEChemicalReaction::FEChemicalReaction(FEModel* pfem) : FEMaterial(pfem)
{
	m_Vovr = false; 
	m_pMP = 0; 
	m_pFwd = 0;
	m_pRev = 0;
}

//-----------------------------------------------------------------------------
void FEChemicalReaction::Init() 
{
	if (m_pFwd) InitializeReactionRate(m_pFwd);
	if (m_pRev) InitializeReactionRate(m_pRev);
}

//-----------------------------------------------------------------------------
void FEChemicalReaction::SetParameter(FEParam& p)
{
	if (strcmp(p.m_szname, "Vbar") == 0)
	{
		m_Vovr = true;
	}
}

//-----------------------------------------------------------------------------
bool FEChemicalReaction::SetParameterAttribute(FEParam& p, const char* szatt, const char* szval)
{
	if (strcmp(p.m_szname, "vR") == 0)
	{
		if (strcmp(szatt, "sbm") == 0)
		{
			int id = atoi(szval) - 1;
			if (id < 0) return false;
			SetSolidReactantsCoefficients(id, m_vRtmp);
		}
		if (strcmp(szatt, "sol") == 0)
		{
			int id = atoi(szval) - 1;
			if ((id < 0) || (id >= MAX_CDOFS)) return false;
			SetSoluteReactantsCoefficients(id, m_vRtmp);
		}
	}
	else if (strcmp(p.m_szname, "vP") == 0)
	{
		if (strcmp(szatt, "sbm") == 0)
		{
			int id = atoi(szval) - 1;
			if (id < 0) return false;
			SetSolidProductsCoefficients(id, m_vPtmp);
		}
		if (strcmp(szatt, "sol") == 0)
		{
			int id = atoi(szval) - 1;
			if ((id < 0) || (id >= MAX_CDOFS)) return false;
			SetSoluteProductsCoefficients(id, m_vPtmp);
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
int FEChemicalReaction::Properties()
{
	return 2;
}

//-----------------------------------------------------------------------------
FEMaterial* FEChemicalReaction::GetProperty(int i)
{
	if (i==0) return m_pFwd;
	if (i==1) return m_pRev;
	return 0;
}

//-----------------------------------------------------------------------------
//! find a material property index ( returns <0 for error)
int FEChemicalReaction::FindPropertyIndex(const char* szname)
{ 
	if (strcmp(szname, "forward_rate") == 0) return 0;
	if (strcmp(szname, "reverse_rate") == 0) return 1;
	return -1; 
}

//-----------------------------------------------------------------------------
//! set a material property (returns false on error)
bool FEChemicalReaction::SetProperty(int i, FEMaterial* pm)
{ 
	if (i==0)
	{
		FEReactionRate* pmr = dynamic_cast<FEReactionRate*>(pm);
		if (pmr) { SetForwardReactionRate(pmr); return true; }
	}
	if (i==1)
	{
		FEReactionRate* pmr = dynamic_cast<FEReactionRate*>(pm);
		if (pmr) { SetReverseReactionRate(pmr); return true; }
	}
	return false; 
}

//-----------------------------------------------------------------------------
// initialize chemical reaction rate
void FEChemicalReaction::InitializeReactionRate(FEReactionRate* m_pRate)
{
	m_pRate->m_pReact = this; 
}
