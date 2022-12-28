/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
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
    int varAC = fedofs.GetVariableIndex("concentration tderiv");
	int cdofs = fedofs.GetVariableSize(varC);
    int ddofs = fedofs.GetVariableSize(varD);
	char sz[8] = {0};
	sprintf(sz, "c%d", cdofs+1);
	fedofs.AddDOF(varC, sz);
    sprintf(sz, "d%d", ddofs+1);
    fedofs.AddDOF(varD, sz);
    sprintf(sz, "ac%d", cdofs+1);
    fedofs.AddDOF(varAC, sz);

	return true;
}

//=============================================================================
// FESolute
//=============================================================================

//-----------------------------------------------------------------------------
// Material parameters for FESoluteData
BEGIN_FECORE_CLASS(FESoluteMaterial, FESolute)

	// These parameters cannot (or should not) be set in the input file since
	// they are copied from the FESoluteData class. 
//	ADD_PARAMETER(m_rhoT, "density");
//	ADD_PARAMETER(m_M, "molar_mass");
//	ADD_PARAMETER(m_z, "charge_number");

	ADD_PARAMETER(m_ID, "sol", FE_PARAM_ATTRIBUTE, "$(solutes)");

	// set material properties
	ADD_PROPERTY(m_pDiff , "diffusivity");
	ADD_PROPERTY(m_pSolub, "solubility");
	ADD_PROPERTY(m_pSupp , "supply", FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESoluteMaterial::FESoluteMaterial(FEModel* fem) : FESolute(fem)
{

}

//-----------------------------------------------------------------------------
//! FESolute constructor

FESolute::FESolute(FEModel* pfem) : FEMaterialProperty(pfem)
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
	if (FEMaterialProperty::Init() == false) return false;

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
	FEMaterialProperty::Serialize(ar);
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
BEGIN_FECORE_CLASS(FESolidBoundMolecule, FEMaterialProperty)

	ADD_PARAMETER(m_ID, "sbm", FE_PARAM_ATTRIBUTE, "$(sbms)");

	ADD_PARAMETER(m_rho0  , "rho0"  );
	ADD_PARAMETER(m_rhomin, "rhomin");
	ADD_PARAMETER(m_rhomax, "rhomax");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! FESolidBoundMolecule constructor

FESolidBoundMolecule::FESolidBoundMolecule(FEModel* pfem) : FEMaterialProperty(pfem)
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
	if (FEMaterialProperty::Init() == false) return false;
	
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
	FEMaterialProperty::Serialize(ar);
	ar & m_ID;
	ar & m_rhoT & m_M & m_z & m_rho0;
	ar & m_rhomin & m_rhomax;
}
