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
    // get number of DOFS
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_CDOFS = fedofs.GetCDOFS();
    
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
int FEChemicalReaction::Properties()
{
	return 2;
}

//-----------------------------------------------------------------------------
FECoreBase* FEChemicalReaction::GetProperty(int i)
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
bool FEChemicalReaction::SetProperty(int i, FECoreBase* pm)
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
    m_pRate->Init();
}

//-----------------------------------------------------------------------------
//! Data serialization
void FEChemicalReaction::Serialize(DumpFile& ar)
{
	FEMaterial::Serialize(ar);
	
	if (ar.IsSaving())
	{
		itrmap p;
		ar << m_nsol << m_vR << m_vP << m_v << m_Vbar << m_Vovr << m_vRtmp << m_vPtmp;
		ar << m_solR.size();
		for (p = m_solR.begin(); p!=m_solR.end(); ++p) {ar << p->first; ar << p->second;}
		ar << m_solP.size();
		for (p = m_solP.begin(); p!=m_solP.end(); ++p) {ar << p->first; ar << p->second;}
		ar << m_sbmR.size();
		for (p = m_sbmR.begin(); p!=m_sbmR.end(); ++p) {ar << p->first; ar << p->second;}
		ar << m_sbmP.size();
		for (p = m_sbmP.begin(); p!=m_sbmP.end(); ++p) {ar << p->first; ar << p->second;}

		if (m_pFwd)
		{
			ar << 1;
			ar << m_pFwd->GetTypeStr(); m_pFwd->Serialize(ar);
		}
		else ar << 0;
		if (m_pRev)
		{
			ar << 1;
			ar << m_pRev->GetTypeStr(); m_pRev->Serialize(ar);
		}
		else ar << 0;

	}
	else
	{
		ar >> m_nsol >> m_vR >> m_vP >> m_v >> m_Vbar >> m_Vovr >> m_vRtmp >> m_vPtmp;
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
		int rr;
		char sz[256] = {0};
		ar >> rr;
		if (rr)
		{
			ar >> sz;
			m_pFwd = dynamic_cast<FEReactionRate*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
			assert(m_pFwd); m_pFwd->Serialize(ar);
			m_pFwd->Init();
		}
		ar >> rr;
		if (rr)
		{
			ar >> sz;
			m_pRev = dynamic_cast<FEReactionRate*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
			assert(m_pRev); m_pRev->Serialize(ar);
			m_pRev->Init();
		}

	}

}