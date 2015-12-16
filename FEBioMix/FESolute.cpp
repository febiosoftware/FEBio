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
    DOFS& fedofs = pfem->GetDOFS();
	int cdofs = fedofs.GetDOFSize("c");
	fedofs.ChangeDOFSize("c", cdofs+1);
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

	// set material properties
	AddProperty(&m_pDiff , "diffusivity");
	AddProperty(&m_pSolub, "solubility" );
	AddProperty(&m_pSupp , "supply"     , false);
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
	
	if (ar.IsSaving())
	{
		ar << GetSoluteID();
	}
	else
	{
		int solID;
		ar >> solID;
		SetSoluteID(solID);
	}
}

//-----------------------------------------------------------------------------
bool FESolute::SetAttribute(const char* szname, const char* szval)
{
    // get number of DOFS
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_CDOFS = fedofs.GetDOFSize("c");
    
	if (strcmp(szname, "sol") == 0)
	{
		int nid = atoi(szval) - 1;
		if ((nid < 0) || (nid >= MAX_CDOFS)) return false;
		SetSoluteID(nid);
	}
	return true;
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
