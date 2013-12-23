#include "FESolutesMaterialPoint.h"

//=============================================================================
//   FESolutesMaterialPoint
//=============================================================================


//-----------------------------------------------------------------------------
//! Create a shallow copy of the material point data
FEMaterialPoint* FESolutesMaterialPoint::Copy()
{
	FESolutesMaterialPoint* pt = new FESolutesMaterialPoint(*this);
	if (m_pt) pt->m_pt = m_pt->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! Initialize material point data
void FESolutesMaterialPoint::Init(bool bflag)
{
	if (bflag)
	{
		m_nsol = m_nsbm = 0;
		m_psi = m_cF = 0;
		m_Ie = vec3d(0,0,0);
		m_rhor = 0;
	}
    
	if (m_pt) m_pt->Init(bflag);
}

//-----------------------------------------------------------------------------
//! Serialize material point data to the archive
void FESolutesMaterialPoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		dmp << m_nsol << m_psi << m_cF << m_Ie << m_nsbm;
		for (int i=0; i<m_nsol; ++i) {
			dmp << m_c[i] << m_gradc[i] << m_j[i] << m_ca[i]
			<< m_k[i] << m_dkdJ[i];
			for (int j=0; j<m_nsol; ++j)
				dmp << m_dkdc[i][j];
		}
		for (int i=0; i<m_nsbm; ++i)
			dmp << m_sbmr[i] << m_sbmrp[i] << m_sbmrhat[i];
	}
	else
	{
		dmp >> m_nsol >> m_psi >> m_cF >> m_Ie >> m_nsbm;
		for (int i=0; i<m_nsol; ++i) {
			dmp >> m_c[i] >> m_gradc[i] >> m_j[i] >> m_ca[i]
			>> m_k[i] >> m_dkdJ[i];
			for (int j=0; j<m_nsol; ++j)
				dmp >> m_dkdc[i][j];
		}
		for (int i=0; i<m_nsbm; ++i)
			dmp >> m_sbmr[i] >> m_sbmrp[i] >> m_sbmrhat[i];
	}
    
	if (m_pt) m_pt->ShallowCopy(dmp, bsave);
}

//-----------------------------------------------------------------------------
//! Serialize material point data to the archive
void FESolutesMaterialPoint::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_nsol << m_psi << m_cF << m_Ie << m_nsbm;
		for (int i=0; i<m_nsol; ++i) {
			ar << m_c[i] << m_gradc[i] << m_j[i] << m_ca[i]
			<< m_k[i] << m_dkdJ[i];
			for (int j=0; j<m_nsol; ++j)
				ar << m_dkdc[i][j];
		}
		for (int i=0; i<m_nsbm; ++i)
			ar << m_sbmr[i] << m_sbmrp[i] << m_sbmrhat[i];
	}
	else
	{
		ar >> m_nsol >> m_psi >> m_cF >> m_Ie >> m_nsbm;
        
		m_c.resize(m_nsol);
		m_gradc.resize(m_nsol);
		m_j.resize(m_nsol);
		m_ca.resize(m_nsol);
		m_k.resize(m_nsol);
		m_dkdJ.resize(m_nsol);
		m_dkdc.resize(m_nsol);
        
		for (int i=0; i<m_nsol; ++i) {
			ar >> m_c[i] >> m_gradc[i] >> m_j[i] >> m_ca[i]
			>> m_k[i] >> m_dkdJ[i];
            
			m_dkdc[i].resize(m_nsol);
			for (int j=0; j<m_nsol; ++j)
			{
				ar >> m_dkdc[i][j];
			}
		}
		
		if (m_nsbm)
		{
			m_sbmr.resize(m_nsbm);
			m_sbmrp.resize(m_nsbm);
			m_sbmrhat.resize(m_nsbm);
            
			for (int i=0; i<m_nsbm; ++i)
				ar >> m_sbmr[i] >> m_sbmrp[i] >> m_sbmrhat[i];
		}
	}
    
	if (m_pt) m_pt->Serialize(ar);
}

