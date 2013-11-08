#include "FESolute.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
//! FESolute constructor

FESolute::FESolute()
{
	m_rhoT = 0;
	m_M = 0;
	m_z = 0;

	AddComponent<FESoluteDiffusivity>(&m_pDiff , "diffusivity");
	AddComponent<FESoluteSolubility >(&m_pSolub, "solubility" );
	AddComponent<FESoluteSupply     >(&m_pSupp , "supply"     );
}

//-----------------------------------------------------------------------------
void FESolute::Init()
{
	FEMaterial::Init();
	m_pDiff->Init();
	m_pSolub->Init();
	if (m_pSupp) m_pSupp->Init();

	FESoluteData* psd = FEModel::FindSD(m_ID);
	if (psd == 0) throw MaterialError("no match with global solute data");
	m_rhoT = psd->m_rhoT;
	m_M = psd->m_M;
	m_z = (int) psd->m_z;
	SetName(psd->m_szname);
	
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
		ar << febio.GetTypeStr<FEMaterial>(m_pSupp ); m_pSupp ->Serialize(ar);
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
		
		ar >> sz;
		m_pSupp = dynamic_cast<FESoluteSupply*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pSupp); m_pSupp->Serialize(ar);
		m_pSupp->Init();
	}
}

//-----------------------------------------------------------------------------
int FESolute::Properties()
{
	return (m_pSupp ? 3 : 2);
}

//-----------------------------------------------------------------------------
FEMaterial* FESolute::GetProperty(int i)
{
	switch (i)
	{
	case 0: return m_pDiff;
	case 1: return m_pSolub;
	case 2: return m_pSupp;
	}
	assert(false);
	return 0;
}

//-----------------------------------------------------------------------------
//! Find the index of a material property
int FESolute::FindPropertyIndex(const char* szname)
{
	if (strcmp(szname, "diffusivity") == 0) return 0;
	if (strcmp(szname, "solubility" ) == 0) return 1;
	if (strcmp(szname, "supply"     ) == 0) return 2;
	return -1;
}

//-----------------------------------------------------------------------------
//! Set a material property
bool FESolute::SetProperty(int n, FEMaterial* pm)
{
	switch(n)
	{
	case 0:
		{
			FESoluteDiffusivity* pmd = dynamic_cast<FESoluteDiffusivity*>(pm);
			if (pmd) { m_pDiff = pmd; return true; }
		}
		break;
	case 1: 
		{
			FESoluteSolubility* pms = dynamic_cast<FESoluteSolubility*>(pm);
			if (pms) { m_pSolub = pms; return true; }
		}
		break;
	case 2:
		{
			FESoluteSupply* pms = dynamic_cast<FESoluteSupply*>(pm);
			if (pms) { m_pSupp = pms; return true; }
		}
		break;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FESolute::SetAttribute(const char* szname, const char* szval)
{
	if (strcmp(szname, "sol") == 0)
	{
		int nid = atoi(szval) - 1;
		if ((nid < 0) || (nid >= MAX_CDOFS)) return false;
		SetSoluteID(nid);
	}
	return true;
}

//-----------------------------------------------------------------------------
FEParam* FESolute::GetParameter(const ParamString& s)
{
	if (s.count() == 1) return FEMultiMaterial::GetParameter(s);
	if      (s == "diffusivity") return m_pDiff ->GetParameter(s.next());
	else if (s == "solubility" ) return m_pSolub->GetParameter(s.next());
	return 0;
}

// Material parameters for the FESolidBoundMolecule material
BEGIN_PARAMETER_LIST(FESolidBoundMolecule, FEMaterial)
	ADD_PARAMETER(m_rho0, FE_PARAM_DOUBLE, "rho0");
	ADD_PARAMETER(m_rhomin, FE_PARAM_DOUBLE, "rhomin");
	ADD_PARAMETER(m_rhomax, FE_PARAM_DOUBLE, "rhomax");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! FESolidBoundMolecule constructor

FESolidBoundMolecule::FESolidBoundMolecule()
{
	m_rhoT = 1;
	m_M = 1;
	m_z = 0;
	m_rho0 = 0;
	m_rhomin = 0;
	m_rhomax = 0;
}

//-----------------------------------------------------------------------------
void FESolidBoundMolecule::Init()
{
	FEMaterial::Init();
	
	FESBMData* psd = FEModel::FindSBM(m_ID);
	if (psd == 0) throw MaterialError("no match with global solid-bound molecule data");
	m_rhoT = psd->m_rhoT;
	m_M = psd->m_M;
	m_z = psd->m_z;
	SetName(psd->m_szname);
	
	if (m_rhoT < 0) throw MaterialError("density must be positive");
	if (m_M < 0) throw MaterialError("molar_mass must be positive");
	
}

//-----------------------------------------------------------------------------
bool FESolidBoundMolecule::SetAttribute(const char* szname, const char* szval)
{
	if (strcmp(szname, "sbm") == 0)
	{
		int nid = atoi(szval) - 1;
		if (nid < 0) return false;
		SetSBMID(nid);
	}
	return true;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FESolidBoundMolecule::Serialize(DumpFile& ar)
{
	FEMaterial::Serialize(ar);
	FEBioKernel& febio = FEBioKernel::GetInstance();
	
	if (ar.IsSaving())
	{
		ar << m_rhoT << m_M << m_z << m_rho0 << m_rhomin << m_rhomax;
	}
	else
	{
		ar >> m_rhoT >> m_M >> m_z >> m_rho0 >> m_rhomin >> m_rhomax;
	}
}

//-----------------------------------------------------------------------------
FEParam* FESolidBoundMolecule::GetParameter(const ParamString& s)
{
	if (s.count() == 1) return FEMaterial::GetParameter(s);
	
	return 0;
}
