#include "stdafx.h"
#include "FEChemicalReaction.h"
#include "FECore/FEElementTraits.h"
#include "FECore/DOFS.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
#include "FEMultiphasic.h"
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
	AddProperty(&m_pFwd ,"forward_rate", 0);
	AddProperty(&m_pRev ,"reverse_rate", 0);
}

//-----------------------------------------------------------------------------
bool FEChemicalReaction::Init() 
{
	// make sure the parent class is set
	assert(m_pMP);
	if (m_pMP == 0) return MaterialError("Parent class not set");

	// set the parents for the reaction rates
	if (m_pFwd) m_pFwd->m_pReact = this;
	if (m_pRev) m_pRev->m_pReact = this;

	// now call base class
	if (FEMaterial::Init() == false) return false;

	// initialize the reaction coefficients
	int isol, isbm, itot;

	const int nsol = m_pMP->Solutes();
	const int nsbm = m_pMP->SBMs();
	const int ntot = nsol + nsbm;

	// initialize the stoichiometric coefficients to zero
	m_nsol = nsol;
	m_vR.assign(ntot, 0);
	m_vP.assign(ntot, 0);
	m_v.assign(ntot, 0);

	// cycle through all the solutes in the mixture and determine
	// if they participate in this reaction
	itrmap it;
	intmap solR = m_solR;
	intmap solP = m_solP;
	for (isol = 0; isol<nsol; ++isol) {
		int sid = m_pMP->GetSolute(isol)->GetSoluteID();
		it = solR.find(sid);
		if (it != solR.end()) m_vR[isol] = it->second;
		it = solP.find(sid);
		if (it != solP.end()) m_vP[isol] = it->second;
	}

	// cycle through all the solid-bound molecules in the mixture
	// and determine if they participate in this reaction
	intmap sbmR = m_sbmR;
	intmap sbmP = m_sbmP;
	for (isbm = 0; isbm<nsbm; ++isbm) {
		int sid = m_pMP->GetSBM(isbm)->GetSBMID();
		it = sbmR.find(sid);
		if (it != sbmR.end()) m_vR[nsol + isbm] = it->second;
		it = sbmP.find(sid);
		if (it != sbmP.end()) m_vP[nsol + isbm] = it->second;
	}

	// evaluate the net stoichiometric coefficient
	for (itot = 0; itot<ntot; ++itot) {
		m_v[itot] = m_vP[itot] - m_vR[itot];
	}

	// evaluate the weighted molar volume of reactants and products
	if (!m_Vovr) {
		m_Vbar = 0;
		for (isol = 0; isol<nsol; ++isol)
			m_Vbar += m_v[isol] * m_pMP->GetSolute(isol)->MolarMass() / m_pMP->GetSolute(isol)->Density();
		for (isbm = 0; isbm<nsbm; ++isbm)
			m_Vbar += m_v[nsol + isbm] * m_pMP->GetSBM(isbm)->MolarMass() / m_pMP->GetSBM(isbm)->Density();
	}

	// check that the chemical reaction satisfies electroneutrality
	int znet = 0;
	for (isol = 0; isol<nsol; ++isol)
		znet += m_v[isol] * m_pMP->GetSolute(isol)->ChargeNumber();
	for (isbm = 0; isbm<nsbm; ++isbm)
		znet += m_v[nsol + isbm] * m_pMP->GetSBM(isbm)->ChargeNumber();
	if (znet != 0) return MaterialError("chemical reaction must satisfy electroneutrality");

	return true;
}

//-----------------------------------------------------------------------------
void FEChemicalReaction::SetParameter(FEParam& p)
{
	if (strcmp(p.name(), "Vbar") == 0)
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
    
	if (strcmp(p.name(), "vR") == 0)
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
	else if (strcmp(p.name(), "vP") == 0)
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

	if (ar.IsShallow() == false)
	{
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
			// restore pointers
			if (m_pFwd) m_pFwd->m_pReact = this;
			if (m_pRev) m_pRev->m_pReact = this;

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
}
