#include "stdafx.h"
#include "FEBioGlobalsSection.h"
#include "FEBioMech/FEElasticMultigeneration.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
//!  This function reads the global variables from the xml file
//!
void FEBioGlobalsSection::Parse(XMLTag& tag)
{
	++tag;
	do
	{
		if      (tag == "Constants"  ) ParseConstants(tag);
		else if (tag == "Solutes"    ) ParseGSSoluteData(tag);
		else if (tag == "SolidBoundMolecules") ParseGSSBMData(tag);
		else if (tag == "Generations") ParseMGData(tag);
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
void FEBioGlobalsSection::ParseGSSoluteData(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	
	// count how many solutes there are
	int nsol = 0;
	XMLTag t(tag); ++t;
	while (!(t == "Solutes")) { 
		if ((t == "solute") && !t.isend()) ++nsol;	// only count on opening tag
		++t;
	}
	
	// read the global solute data
	++tag;
	for (int i=0; i<nsol; ++i)
	{
		FESoluteData* psd = new FESoluteData;
		psd->m_nID = atoi(tag.AttributeValue("id"))-1;
		const char* sz = tag.AttributeValue("name");
		
		if (strcmp(sz, "") == 0)
			throw XMLReader::InvalidAttributeValue(tag, "name", sz);
		
		strcpy(psd->m_szname, sz);

		// read solute properties
		++tag;
		do
		{
/*			if (tag == "charge_number")
			{
				tag.value(psd->m_z);
			}
			else if (tag == "density")
			{
				tag.value(psd->m_rhoT);
			}
			else if (tag == "molar_mass")
			{
				tag.value(psd->m_M);
			}
			else throw XMLReader::InvalidTag(tag);
*/			
			if (m_pim->ReadParameter(tag, psd->GetParameterList()) == false)
			{
				throw XMLReader::InvalidTag(tag);
			}
			++tag;
		}
		while (!tag.isend());
		
		fem.AddSoluteData(psd);
		
		++tag;
	}
}

//-----------------------------------------------------------------------------
void FEBioGlobalsSection::ParseGSSBMData(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();

	// count how many solid-bound molecules there are
	int nsbm = 0;
	XMLTag t(tag); ++t;
	while (!(t == "SolidBoundMolecules")) { 
		if ((t == "solid_bound") && !t.isend()) ++nsbm;	// only count on opening tag
		++t;
	}
	
	// read the global sbm data
	++tag;
	for (int i=0; i<nsbm; ++i)
	{
		FESBMData* psd = new FESBMData;
		psd->m_nID = atoi(tag.AttributeValue("id"))-1;
		const char* sz = tag.AttributeValue("name");
		
		if (strcmp(sz, "") == 0)
			throw XMLReader::InvalidAttributeValue(tag, "name", sz);
		
		strcpy(psd->m_szname, sz);
		
		// read sbm properties
		++tag;
		do
		{
			if (tag == "charge_number")
			{
				tag.value(psd->m_z);
			}
			else if (tag == "density")
			{
				tag.value(psd->m_rhoT);
			}
			else if (tag == "molar_mass")
			{
				tag.value(psd->m_M);
			}
			else throw XMLReader::InvalidTag(tag);
			
			++tag;
		}
		while (!tag.isend());
		
		fem.AddSBMData(psd);
		
		++tag;
	}
}

//-----------------------------------------------------------------------------
//! Parse the time increments for multigeneration materials
void FEBioGlobalsSection::ParseMGData(XMLTag &tag)
{
	++tag;
	do
	{
		if (tag == "gen") {
			FEGenerationData* G = new FEGenerationData;
			int id = atoi(tag.AttributeValue("id"))-1;
			if (id) {
				G->born = false;
			} else {
				G->born = true;
			}
			tag.value(G->btime);
			FEElasticMultigeneration::PushGeneration(G);
		}
		++tag;
	}
	while (!tag.isend());

}
