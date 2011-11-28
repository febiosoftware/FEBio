#include "stdafx.h"
#include "FESolute.h"
#include "fem.h"

// register the material with the framework
REGISTER_MATERIAL(FESolute, "solute");

// Material parameters for the FESolute material
BEGIN_PARAMETER_LIST(FESolute, FEMaterial)
ADD_PARAMETER(m_rhoT, FE_PARAM_DOUBLE, "density");
ADD_PARAMETER(m_M, FE_PARAM_DOUBLE, "molar_mass");
ADD_PARAMETER(m_z, FE_PARAM_INT, "charge_number");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! FESolute constructor

FESolute::FESolute()
{	m_pDiff = 0;
	m_pSolub = 0;
	m_rhoT = 0;
	m_M = 0;
	m_z = 0;
}

//-----------------------------------------------------------------------------
void FESolute::Init()
{
	FEMaterial::Init();
	m_pDiff->Init();
	m_pSolub->Init();
	
	if (m_rhoT < 0) throw MaterialError("density must be positive");
	if (m_M < 0) throw MaterialError("molar_mass must be positive");
		
}

//-----------------------------------------------------------------------------
//! Data serialization
void FESolute::Serialize(DumpFile& ar)
{
	FEMaterial::Serialize(ar);
	FEBioKernel& febio = FEBioKernel::GetInstance();
	
	if (ar.IsSaving())
	{
		ar << m_rhoT << m_M << m_z;
		
		ar << febio.GetTypeStr<FEMaterial>(m_pDiff ); m_pDiff ->Serialize(ar);
		ar << febio.GetTypeStr<FEMaterial>(m_pSolub); m_pSolub->Serialize(ar);
	}
	else
	{
		ar >> m_rhoT >> m_M >> m_z;
		
		char sz[256] = {0};
		ar >> sz;
		m_pDiff = dynamic_cast<FESoluteDiffusivity*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pDiff); m_pDiff->Serialize(ar);
		m_pDiff->Init();
		
		ar >> sz;
		m_pSolub = dynamic_cast<FESoluteSolubility*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pSolub); m_pSolub->Serialize(ar);
		m_pSolub->Init();
		
	}
}
