#include "stdafx.h"
#include "FESolutesMaterialPoint.h"
#include "FECore/DumpStream.h"

//=============================================================================
//   FESolutesMaterialPoint
//=============================================================================


//-----------------------------------------------------------------------------
//! Create a shallow copy of the material point data
FEMaterialPoint* FESolutesMaterialPoint::Copy()
{
	FESolutesMaterialPoint* pt = new FESolutesMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! Initialize material point data
void FESolutesMaterialPoint::Init()
{
	m_nsol = m_nsbm = 0;
	m_psi = m_cF = 0;
	m_Ie = vec3d(0,0,0);
	m_rhor = 0;
    m_c.clear();
    m_gradc.clear();
    m_j.clear();
    m_ca.clear();
    m_crp.clear();
    m_sbmr.clear();
    m_sbmrp.clear();
    m_sbmrhat.clear();
    m_sbmrhatp.clear();
    m_sbmrmin.clear();
    m_sbmrmax.clear();
    m_k.clear();
    m_dkdJ.clear();
    m_dkdJJ.clear();
    m_dkdc.clear();
    m_dkdJc.clear();
    m_dkdcc.clear();
    m_dkdr.clear();
    m_dkdJr.clear();
    m_dkdrc.clear();
    m_cri.clear();
    m_crd.clear();
    m_strain = 0;
    m_pe = m_pi = 0;
    m_ce.clear();
    m_ci.clear();
    m_ide.clear();
    m_idi.clear();
    
	// don't forget to initialize the base class
    FEMaterialPoint::Init();
}

//-----------------------------------------------------------------------------
//! Serialize material point data to the archive
void FESolutesMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPoint::Serialize(ar);
	ar & m_nsol & m_psi & m_cF & m_Ie & m_nsbm;
	ar & m_c & m_gradc & m_j & m_ca & m_crp & m_k & m_dkdJ;
	ar & m_dkdc;
	ar & m_sbmr & m_sbmrp & m_sbmrhat & m_sbmrhatp;
	ar & m_cri;
	ar & m_crd;
	ar & m_strain & m_pe & m_pi;
	ar & m_ce & m_ide;
	ar & m_ci & m_idi;
}
