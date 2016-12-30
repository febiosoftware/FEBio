#include "stdafx.h"
#include "FEBioGlobalsSection.h"
#include "FECore/FEModel.h"
#include "FECore/FEGlobalData.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
//!  This function reads the global variables from the xml file
//!
void FEBioGlobalsSection::Parse(XMLTag& tag)
{
	++tag;
	do
	{
		if      (tag == "Constants"          ) ParseConstants(tag);
		else if (tag == "Solutes"            ) ParseGlobalData(tag);
		else if (tag == "SolidBoundMolecules") ParseGlobalData(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioGlobalsSection::ParseConstants(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	++tag;
	string s;
	double v;
	do
	{
		s = string(tag.Name());
		tag.value(v);
		fem.SetGlobalConstant(s, v);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioGlobalsSection::ParseGlobalData(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	
	// read the global solute data
	++tag;
	do
	{
		// create new global data
		FEGlobalData* pgd = fecore_new<FEGlobalData>(FEGLOBALDATA_ID, tag.Name(), &fem);

		// TODO: We have to call the Init member here because solute data 
		//       allocates the concentration dofs and they have to be allocated before 
		//       materials are read in. I'd like to move this to FEModel::Init but not sure
		//       yet how. 
		pgd->Init();

		// assign attributes
		int natt = tag.m_natt;
		for (int i=0; i<natt; ++i) pgd->SetAttribute(tag.m_att[i].m_szatt, tag.m_att[i].m_szatv);

		// read solute properties
		if (tag.isleaf() == false)
		{
			FEParameterList& pl = pgd->GetParameterList();
			++tag;
			do
			{
				if (m_pim->ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
		
		fem.AddGlobalData(pgd);
		
		++tag;
	}
	while (!tag.isend());
}
