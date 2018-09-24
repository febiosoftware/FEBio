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
	ADD_PARAMETER(m_rhoT, "density");
	ADD_PARAMETER(m_M, "molar_mass");
	ADD_PARAMETER(m_z, "charge_number");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FESoluteData::FESoluteData(FEModel* pfem) : FEGlobalData(pfem)
{ 
	m_rhoT = 1; 
	m_M = 1; 
	m_z = 0; 
}

//-----------------------------------------------------------------------------
// TODO: Maybe I can use the ID to make sure the dof is not duplicated.
bool FESoluteData::Init()
{
	// for each solute we have to add a concentration degree of freedom
	FEModel& fem = *GetFEModel();
    DOFS& fedofs = fem.GetDOFS();
	int varC = fedofs.GetVariableIndex("concentration");
    int varD = fedofs.GetVariableIndex("shell concentration");
	int cdofs = fedofs.GetVariableSize(varC);
    int ddofs = fedofs.GetVariableSize(varD);
	char sz[8] = {0};
	sprintf(sz, "c%d", cdofs+1);
	fedofs.AddDOF(varC, sz);
    sprintf(sz, "d%d", ddofs+1);
    fedofs.AddDOF(varD, sz);

	return true;
}

//-----------------------------------------------------------------------------
bool FESoluteData::SetAttribute(const char* szname, const char* szval)
{
	if (strcmp(szname, "id") == 0)
	{
		SetID(atoi(szval)-1);
		return true;
	}
	else if (strcmp(szname, "name") == 0)
	{
		if (strcmp(szval, "") == 0) return false;
		SetName(szval);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Store the solute data to the archive
void FESoluteData::Serialize(DumpStream &ar)
{
	if (ar.IsSaving()) ar << GetID();
	else 
	{
		int nid;
		ar >> nid;
		SetID(nid);
	}

	// store parameters
	FEGlobalData::Serialize(ar);
}

//=============================================================================
// FESolute
//=============================================================================

//-----------------------------------------------------------------------------
// Material parameters for FESoluteData
BEGIN_PARAMETER_LIST(FESolute, FEMaterial)
	ADD_PARAMETER(m_rhoT, "density");
	ADD_PARAMETER(m_M, "molar_mass");
	ADD_PARAMETER(m_z, "charge_number");
END_PARAMETER_LIST();


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
	AddProperty(&m_pSupp , "supply"     , 0);
}

//-----------------------------------------------------------------------------
FESoluteData* FESolute::FindSoluteData(int nid)
{
	FEModel& fem = *GetFEModel();
	int N = GetFEModel()->GlobalDataItems();
	for (int i=0; i<N; ++i)
	{
		FESoluteData* psd = dynamic_cast<FESoluteData*>(fem.GetGlobalData(i));
		if (psd && (psd->GetID() == nid)) return psd;
	}
	return 0;
}

//-----------------------------------------------------------------------------
bool FESolute::Init()
{
	if (FEMaterial::Init() == false) return false;

	FESoluteData* psd = FindSoluteData(m_ID);
	if (psd == 0) return MaterialError("no match with global solute data");
	m_rhoT = psd->m_rhoT;
	m_M = psd->m_M;
	m_z = (int) psd->m_z;
	SetName(psd->GetName());
	
	if (m_rhoT < 0) return MaterialError("density must be positive");
	if (m_M < 0) return MaterialError("molar_mass must be positive");		

	return true;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FESolute::Serialize(DumpStream& ar)
{
	FEMaterial::Serialize(ar);
	
	if (ar.IsSaving())
	{
		ar << GetSoluteID();
		ar << GetSoluteLocalID();
	}
	else
	{
		int solID, solLID;
		ar >> solID >> solLID;
		SetSoluteID(solID);
		SetSoluteLocalID(solLID);
	}
}

//-----------------------------------------------------------------------------
bool FESolute::SetAttribute(const char* szname, const char* szval)
{
    // get number of DOFS
    DOFS& fedofs = GetFEModel()->GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize("concentration");
    
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
	ADD_PARAMETER(m_rhoT, "density"      );
	ADD_PARAMETER(m_M   , "molar_mass"   );
	ADD_PARAMETER(m_z   , "charge_number");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FESBMData::FESBMData(FEModel* pfem) : FEGlobalData(pfem)
{ 
	m_rhoT = 1; 
	m_M = 1; 
	m_z = 0; 
}

//-----------------------------------------------------------------------------
bool FESBMData::SetAttribute(const char* szname, const char* szval)
{
	if (strcmp(szname, "id") == 0)
	{
		SetID(atoi(szval)-1);
		return true;
	}
	else if (strcmp(szname, "name") == 0)
	{
		if (strcmp(szval, "") == 0) return false;
		SetName(szval);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Store the solute data to the archive
void FESBMData::Serialize(DumpStream& ar)
{
	if (ar.IsSaving()) ar << GetID();
	else 
	{
		int nid;
		ar >> nid;
		SetID(nid);
	}

	// store parameters
	FEGlobalData::Serialize(ar);
}

//=============================================================================
// FESolidBoundMolecule
//=============================================================================

// Material parameters for the FESolidBoundMolecule material
BEGIN_PARAMETER_LIST(FESolidBoundMolecule, FEMaterial)
	ADD_PARAMETER(m_rho0  , "rho0"  );
	ADD_PARAMETER(m_rhomin, "rhomin");
	ADD_PARAMETER(m_rhomax, "rhomax");
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
		if (psd && (psd->GetID() == nid)) return psd;
	}
	return 0;
}

//-----------------------------------------------------------------------------
bool FESolidBoundMolecule::Init()
{
	if (FEMaterial::Init() == false) return false;
	
	FESBMData* psd = FindSBMData(m_ID);
	if (psd == 0) return MaterialError("no match with global solid-bound molecule data");
	m_rhoT = psd->m_rhoT;
	m_M = psd->m_M;
	m_z = psd->m_z;
	SetName(psd->GetName());
	
	if (m_rhoT < 0) return MaterialError("density must be positive");
	if (m_M < 0) return MaterialError("molar_mass must be positive");
	
	return true;
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
void FESolidBoundMolecule::Serialize(DumpStream& ar)
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
