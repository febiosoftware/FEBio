#include "FESolute.h"
#include <FECore/FEModel.h>
#include <FECore/DOFS.h>
#include <FECore/log.h>

//=============================================================================
// FESoluteData
//=============================================================================

//-----------------------------------------------------------------------------
// Material parameters for FESoluteData
BEGIN_FECORE_CLASS(FESoluteData, FEGlobalData)
	ADD_PARAMETER(m_rhoT, "density");
	ADD_PARAMETER(m_M, "molar_mass");
	ADD_PARAMETER(m_z, "charge_number");
END_FECORE_CLASS();

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

//=============================================================================
// FESolute
//=============================================================================

//-----------------------------------------------------------------------------
// Material parameters for FESoluteData
BEGIN_FECORE_CLASS(FESolute, FEMaterial)
	ADD_PARAMETER(m_rhoT, "density");
	ADD_PARAMETER(m_M, "molar_mass");
	ADD_PARAMETER(m_z, "charge_number");

	ADD_PARAMETER(m_ID, "sol", FE_PARAM_ATTRIBUTE, 0);

	// set material properties
	ADD_PROPERTY(m_pDiff , "diffusivity");
	ADD_PROPERTY(m_pSolub, "solubility");
	ADD_PROPERTY(m_pSupp , "supply", FEProperty::Optional);

END_FECORE_CLASS();


//-----------------------------------------------------------------------------
//! FESolute constructor

FESolute::FESolute(FEModel* pfem) : FEMaterial(pfem)
{
	m_rhoT = 0;
	m_M = 0;
	m_z = 0;

	m_ID = -1;

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
		if (psd && (psd->GetID() == nid)) return psd;
	}
	return 0;
}

//-----------------------------------------------------------------------------
bool FESolute::Init()
{
	if (FEMaterial::Init() == false) return false;

	FESoluteData* psd = FindSoluteData(m_ID);
	if (psd == 0) {
		feLogError("no match with global solute data");
		return false;
	}
	m_rhoT = psd->m_rhoT;
	m_M = psd->m_M;
	m_z = (int) psd->m_z;
	SetName(psd->GetName());
	
	if (m_rhoT < 0) { feLogError("density must be positive"   ); return false; }
	if (m_M    < 0) { feLogError("molar_mass must be positive"); return false; }

	return true;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FESolute::Serialize(DumpStream& ar)
{
	FEMaterial::Serialize(ar);
	ar & m_ID & m_LID;
}

//=============================================================================
// FESBMData
//=============================================================================

//-----------------------------------------------------------------------------
// Material parameters for FESoluteData
BEGIN_FECORE_CLASS(FESBMData, FEGlobalData)
	ADD_PARAMETER(m_rhoT, "density"      );
	ADD_PARAMETER(m_M   , "molar_mass"   );
	ADD_PARAMETER(m_z   , "charge_number");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESBMData::FESBMData(FEModel* pfem) : FEGlobalData(pfem)
{ 
	m_rhoT = 1; 
	m_M = 1; 
	m_z = 0; 
}

//=============================================================================
// FESolidBoundMolecule
//=============================================================================

// Material parameters for the FESolidBoundMolecule material
BEGIN_FECORE_CLASS(FESolidBoundMolecule, FEMaterial)
	ADD_PARAMETER(m_rho0  , "rho0"  );
	ADD_PARAMETER(m_rhomin, "rhomin");
	ADD_PARAMETER(m_rhomax, "rhomax");

	ADD_PARAMETER(m_ID, "sbm", FE_PARAM_ATTRIBUTE, 0);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! FESolidBoundMolecule constructor

FESolidBoundMolecule::FESolidBoundMolecule(FEModel* pfem) : FEMaterial(pfem)
{
	m_ID = -1;
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
	if (psd == 0) {
		feLogError("no match with global solid-bound molecule data");
		return false;
	}
	m_rhoT = psd->m_rhoT;
	m_M = psd->m_M;
	m_z = psd->m_z;
	SetName(psd->GetName());
	
	if (m_rhoT < 0) { feLogError("density must be positive"   ); return false; }
	if (m_M    < 0) { feLogError("molar_mass must be positive"); return false; }
	
	return true;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FESolidBoundMolecule::Serialize(DumpStream& ar)
{
	FEMaterial::Serialize(ar);
	ar & m_ID;
	ar & m_rhoT & m_M & m_z & m_rho0;
	ar & m_rhomin & m_rhomax;
}
