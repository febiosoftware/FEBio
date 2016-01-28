#include "stdafx.h"
#include "FEChemicalReaction.h"
#include "FECore/FEElementTraits.h"
#include "FECore/DOFS.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
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

	// set material properties
	AddProperty(&m_pFwd ,"forward_rate", false);
	AddProperty(&m_pRev ,"reverse_rate", false);
}

//-----------------------------------------------------------------------------
bool FEChemicalReaction::Init() 
{
	if (m_pFwd) m_pFwd->m_pReact = this;
	if (m_pRev) m_pRev->m_pReact = this;
	return FEMaterial::Init();
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
    // get number of DOFS
    DOFS& fedofs = GetFEModel()->GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize("concentration");
    
	if (strcmp(p.m_szname, "vR") == 0)
	{
		if (strcmp(szatt, "sbm") == 0)
		{
			int id = atoi(szval) - 1;
			if (id < 0) return false;
			SetSolidReactantsCoefficients(id, m_vRtmp);
			return true;
		}
		if (strcmp(szatt, "sol") == 0)
		{
			int id = atoi(szval) - 1;
			if ((id < 0) || (id >= MAX_CDOFS)) return false;
			SetSoluteReactantsCoefficients(id, m_vRtmp);
			return true;
		}
	}
	else if (strcmp(p.m_szname, "vP") == 0)
	{
		if (strcmp(szatt, "sbm") == 0)
		{
			int id = atoi(szval) - 1;
			if (id < 0) return false;
			SetSolidProductsCoefficients(id, m_vPtmp);
			return true;
		}
		if (strcmp(szatt, "sol") == 0)
		{
			int id = atoi(szval) - 1;
			if ((id < 0) || (id >= MAX_CDOFS)) return false;
			SetSoluteProductsCoefficients(id, m_vPtmp);
			return true;
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FEChemicalReaction::Serialize(DumpStream& ar)
{
	FEMaterial::Serialize(ar);
	
	if (ar.IsSaving())
	{
		itrmap p;
		ar << m_nsol << m_vR << m_vP << m_v << m_Vovr;
		ar << (int) m_solR.size();
		for (p = m_solR.begin(); p!=m_solR.end(); ++p) {ar << p->first; ar << p->second;}
		ar << (int) m_solP.size();
		for (p = m_solP.begin(); p!=m_solP.end(); ++p) {ar << p->first; ar << p->second;}
		ar << (int) m_sbmR.size();
		for (p = m_sbmR.begin(); p!=m_sbmR.end(); ++p) {ar << p->first; ar << p->second;}
		ar << (int) m_sbmP.size();
		for (p = m_sbmP.begin(); p!=m_sbmP.end(); ++p) {ar << p->first; ar << p->second;}
	}
	else
	{
		ar >> m_nsol >> m_vR >> m_vP >> m_v >> m_Vovr;
		int size, id, vR;
		ar >> size;
		for (int i=0; i<size; ++i)
		{
			ar >> id; ar >> vR;
			SetSoluteReactantsCoefficients(id, vR);
		}
		ar >> size;
		for (int i=0; i<size; ++i)
		{
			ar >> id; ar >> vR;
			SetSoluteProductsCoefficients(id, vR);
		}
		ar >> size;
		for (int i=0; i<size; ++i)
		{
			ar >> id; ar >> vR;
			SetSolidReactantsCoefficients(id, vR);
		}
		ar >> size;
		for (int i=0; i<size; ++i)
		{
			ar >> id; ar >> vR;
			SetSolidProductsCoefficients(id, vR);
		}
	}
}
