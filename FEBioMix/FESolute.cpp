#include "FESolute.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
#include "FECore/DOFS.h"

//=============================================================================
// FESoluteData
//=============================================================================

//-----------------------------------------------------------------------------
// Material parameters for FESoluteData
BEGIN_PARAMETER_LIST(FESoluteData, FEGlobalData)
	ADD_PARAMETER(m_rhoT, FE_PARAM_DOUBLE, "density");
	ADD_PARAMETER(m_M, FE_PARAM_DOUBLE, "molar_mass");
	ADD_PARAMETER(m_z, FE_PARAM_INT, "charge_number");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FESoluteData::FESoluteData(FEModel* pfem)
{ 
	m_nID = -1; 
	m_rhoT = 1; 
	m_M = 1; 
	m_z = 0; 
	m_szname[0] = 0; 

	// for each solute we have to increase the number of concentration dofs
    DOFS& fedofs = *DOFS::GetInstance();
	int ndofs = fedofs.GetNDOFS();
	int cdofs = fedofs.GetCDOFS();
    fedofs.SetCDOFS(cdofs + 1);
    fedofs.SetNDOFS(ndofs + 1);
}

//-----------------------------------------------------------------------------
bool FESoluteData::SetAttribute(const char* szname, const char* szval)
{
	if (strcmp(szname, "id") == 0)
	{
		m_nID = atoi(szval)-1;
		return true;
	}
	else if (strcmp(szname, "name") == 0)
	{
		if (strcmp(szval, "") == 0) return false;
		strcpy(m_szname, szval);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Store the solute data to the archive
void FESoluteData::Serialize(DumpFile &ar)
{
	if (ar.IsSaving()) ar << m_nID;
	else ar >> m_nID;

	// store parameters
	FEGlobalData::Serialize(ar);
}

//=============================================================================
// FESolute
//=============================================================================

//-----------------------------------------------------------------------------
//! FESolute constructor

FESolute::FESolute(FEModel* pfem) : FEMaterial(pfem)
{
	m_rhoT = 0;
	m_M = 0;
	m_z = 0;

	m_pDiff = 0;
	m_pSolub = 0;
	m_pSupp = 0;
}


//-----------------------------------------------------------------------------
FESoluteData* FESolute::FindSoluteData(int nid)
{
	FEModel& fem = *GetFEModel();
	int N = GetFEModel()->GlobalDataItems();
	for (int i=0; i<N; ++i)
	{
		FESoluteData* psd = dynamic_cast<FESoluteData*>(fem.GetGlobalData(i));
		if (psd && (psd->m_nID == nid)) return psd;
	}
	return 0;
}

//-----------------------------------------------------------------------------
void FESolute::Init()
{
	FEMaterial::Init();
	m_pDiff->SetParent(this); m_pDiff->Init();
	m_pSolub->SetParent(this); m_pSolub->Init();
	if (m_pSupp) { m_pSupp->SetParent(this); m_pSupp->Init(); }

	FESoluteData* psd = FindSoluteData(m_ID);
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
	
	int nSupp = 0;
	if (ar.IsSaving())
	{
		ar << GetSoluteID();
		ar << m_rhoT << m_M << m_z;
		
		ar << m_pDiff->GetTypeStr(); m_pDiff ->Serialize(ar);
		ar << m_pSolub->GetTypeStr(); m_pSolub->Serialize(ar);
		
		if (m_pSupp == 0) ar << nSupp;
		else
		{
			nSupp = 1;
			ar << nSupp;
			ar << m_pSupp->GetTypeStr();
			m_pSupp ->Serialize(ar);
		}
	}
	else
	{
		int solID;
		ar >> solID;
		SetSoluteID(solID);
		ar >> m_rhoT >> m_M >> m_z;
		
		char sz[256] = {0};
		ar >> sz;
		m_pDiff = dynamic_cast<FESoluteDiffusivity*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
		assert(m_pDiff); m_pDiff->Serialize(ar);
		m_pDiff->Init();
		
		ar >> sz;
		m_pSolub = dynamic_cast<FESoluteSolubility*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
		assert(m_pSolub); m_pSolub->Serialize(ar);
		m_pSolub->Init();
		
		ar >> nSupp;
		if (nSupp)
		{
			ar >> sz;
			m_pSupp = dynamic_cast<FESoluteSupply*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
			assert(m_pSupp); m_pSupp->Serialize(ar);
			m_pSupp->Init();
		}
	}
}

//-----------------------------------------------------------------------------
int FESolute::Properties()
{
	return (m_pSupp ? 3 : 2);
}

//-----------------------------------------------------------------------------
FECoreBase* FESolute::GetProperty(int i)
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
bool FESolute::SetProperty(int n, FECoreBase* pm)
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
    // get number of DOFS
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_CDOFS = fedofs.GetCDOFS();
    
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
	if (s.count() == 1) return FEMaterial::GetParameter(s);
	if      (s == "diffusivity") return m_pDiff ->GetParameter(s.next());
	else if (s == "solubility" ) return m_pSolub->GetParameter(s.next());
	return 0;
}

//=============================================================================
// FESBMData
//=============================================================================

//-----------------------------------------------------------------------------
// Material parameters for FESoluteData
BEGIN_PARAMETER_LIST(FESBMData, FEGlobalData)
	ADD_PARAMETER(m_rhoT, FE_PARAM_DOUBLE, "density"      );
	ADD_PARAMETER(m_M   , FE_PARAM_DOUBLE, "molar_mass"   );
	ADD_PARAMETER(m_z   , FE_PARAM_INT   , "charge_number");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FESBMData::FESBMData(FEModel* pfem)
{ 
	m_nID = -1; 
	m_rhoT = 1; 
	m_M = 1; 
	m_z = 0; 
	m_szname[0] = 0; 
}

//-----------------------------------------------------------------------------
bool FESBMData::SetAttribute(const char* szname, const char* szval)
{
	if (strcmp(szname, "id") == 0)
	{
		m_nID = atoi(szval)-1;
		return true;
	}
	else if (strcmp(szname, "name") == 0)
	{
		if (strcmp(szval, "") == 0) return false;
		strcpy(m_szname, szval);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Store the solute data to the archive
void FESBMData::Serialize(DumpFile &ar)
{
	if (ar.IsSaving()) ar << m_nID;
	else ar >> m_nID;

	// store parameters
	FEGlobalData::Serialize(ar);
}

//=============================================================================
// FESolidBoundMolecule
//=============================================================================

// Material parameters for the FESolidBoundMolecule material
BEGIN_PARAMETER_LIST(FESolidBoundMolecule, FEMaterial)
	ADD_PARAMETER(m_rho0  , FE_PARAM_DOUBLE, "rho0"  );
	ADD_PARAMETER(m_rhomin, FE_PARAM_DOUBLE, "rhomin");
	ADD_PARAMETER(m_rhomax, FE_PARAM_DOUBLE, "rhomax");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! FESolidBoundMolecule constructor

FESolidBoundMolecule::FESolidBoundMolecule(FEModel* pfem) : FEMaterial(pfem)
{
	m_rhoT = 1;
	m_M = 1;
	m_z = 0;
	m_rho0 = 0;
	m_rhomin = 0;
	m_rhomax = 0;
}

//-----------------------------------------------------------------------------
FESBMData* FESolidBoundMolecule::FindSBMData(int nid)
{
	FEModel& fem = *GetFEModel();
	int N = GetFEModel()->GlobalDataItems();
	for (int i=0; i<N; ++i)
	{
		FESBMData* psd = dynamic_cast<FESBMData*>(fem.GetGlobalData(i));
		if (psd && (psd->m_nID == nid)) return psd;
	}
	return 0;
}

//-----------------------------------------------------------------------------
void FESolidBoundMolecule::Init()
{
	FEMaterial::Init();
	
	FESBMData* psd = FindSBMData(m_ID);
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
	FECoreKernel& febio = FECoreKernel::GetInstance();
	
	if (ar.IsSaving())
	{
		ar << GetSBMID();
		ar << m_rhoT << m_M << m_z << m_rho0 << m_rhomin << m_rhomax;
	}
	else
	{
		int SBMID;
		ar >> SBMID;
		SetSBMID(SBMID);
		ar >> m_rhoT >> m_M >> m_z >> m_rho0 >> m_rhomin >> m_rhomax;
	}
}

//-----------------------------------------------------------------------------
FEParam* FESolidBoundMolecule::GetParameter(const ParamString& s)
{
	if (s.count() == 1) return FEMaterial::GetParameter(s);
	
	return 0;
}
