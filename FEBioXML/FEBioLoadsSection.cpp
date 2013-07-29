#include "stdafx.h"
#include "FEBioLoadsSection.h"
#include <FEBioMech/FEPointBodyForce.h>
#include <FEBioHeat/FEHeatSource.h>

//-----------------------------------------------------------------------------
//!  Parses the loads section from the xml file (version 1.2 or up)
//!
void FEBioLoadsSection::Parse(XMLTag& tag)
{
	assert(m_pim->Version() >= 0x0102);
	
	// make sure this tag has children
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "force"              ) ParseBCForce             (tag);
		else if (tag == "pressure"           ) ParseBCPressure          (tag);
		else if (tag == "traction"           ) ParseBCTraction          (tag);
		else if (tag == "normal_traction"    ) ParseBCPoroNormalTraction(tag);
		else if (tag == "fluidflux"          ) ParseBCFluidFlux         (tag);
		else if (tag == "soluteflux"         ) ParseBCSoluteFlux        (tag);
		else if (tag == "heatflux"           ) ParseBCHeatFlux          (tag);
		else if (tag == "convective_heatflux") ParseBCConvectiveHeatFlux(tag);
		else if (tag == "body_force"         ) ParseBodyForce           (tag);
		else if (tag == "heat_source"        ) ParseHeatSource          (tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
// NOTE: note that this section used to be in the Globals section (version 1.1)
void FEBioLoadsSection::ParseBodyForce(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();

	const char* szt = tag.AttributeValue("type", true);
	if (szt == 0) szt = "const";

	if (strcmp(szt, "point") == 0)
	{
		FEPointBodyForce* pf = new FEPointBodyForce(&fem);
		FEParameterList& pl = pf->GetParameterList();
		++tag;
		do
		{
			if (tag == "a")
			{
				const char* szlc = tag.AttributeValue("lc");
//						pf->lc[0] = pf->lc[1] = pf->lc[2] = atoi(szlc);
				tag.value(pf->m_a);
			}
			else if (tag == "node")
			{
				tag.value(pf->m_inode); 
				pf->m_inode -= 1;
			}
			else if (m_pim->ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());

		fem.AddBodyLoad(pf);
	}
	else
	{
		// see if the kernel knows this force
		FEBioKernel& febio = FEBioKernel::GetInstance();
		FEBodyForce* pf = febio.Create<FEBodyForce>(szt, &fem);
		if (pf)
		{
			if (!tag.isleaf())
			{
				FEParameterList& pl = pf->GetParameterList();
				++tag;
				do
				{
					if (m_pim->ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
					++tag;
				}
				while (!tag.isend());
			}

			fem.AddBodyLoad(pf);
		}
		else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseHeatSource(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEHeatSource* phs = new  FEHeatSource(&fem);
	FEParameterList& PL = phs->GetParameterList();
	++tag;
	do
	{
		if (m_pim->ReadParameter(tag, PL) == false) throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
	fem.AddBodyLoad(phs);
}
