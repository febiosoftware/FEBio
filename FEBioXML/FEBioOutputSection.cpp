#include "stdafx.h"
#include "FEBioOutputSection.h"
#include "FECore/NodeDataRecord.h"
#include "FECore/ElementDataRecord.h"
#include "FECore/ObjectDataRecord.h"
#include "FEBioPlot/FEBioPlotFile.h"
#include "FEBioPlot/LSDYNAPlotFile.h"

//-----------------------------------------------------------------------------
void FEBioOutputSection::Parse(XMLTag& tag)
{
	++tag;
	do
	{
		if (tag == "logfile") ParseLogfile(tag);
		else if (tag == "plotfile") ParsePlotfile(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioOutputSection::ParseLogfile(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	const char* sz;

	// see if the log file has any attributes
	const char* szlog = tag.AttributeValue("file", true);
	if (szlog) m_pim->SetLogfileName(szlog);

	if (tag.isleaf()) return;
	++tag;
	do
	{
		if (tag == "node_data")
		{
			sz = tag.AttributeValue("file", true);

			NodeDataRecord* prec = new NodeDataRecord(&fem, sz);
			const char* szdata = tag.AttributeValue("data");
			prec->Parse(szdata);

			const char* szname = tag.AttributeValue("name", true);
			if (szname != 0) prec->SetName(szname); else prec->SetName(szdata);

			sz = tag.AttributeValue("delim", true);
			if (sz != 0) prec->SetDelim(sz);

			sz = tag.AttributeValue("comments", true);
			if (sz != 0)
			{
				if      (strcmp(sz, "on") == 0) prec->SetComments(true);
				else if (strcmp(sz, "off") == 0) prec->SetComments(false); 
			}

			if (tag.isleaf()) prec->DataRecord::SetItemList(tag.szvalue());
			else
			{
				++tag;
				if (tag == "node_set")
				{
					FENodeSet* pns = 0;
					const char* szid = tag.AttributeValue("id", true);
					if (szid == 0)
					{
						const char* szname = tag.AttributeValue("name");
						pns = mesh.FindNodeSet(szname);
					}
					else pns = mesh.FindNodeSet(atoi(szid));

					if (pns == 0) throw XMLReader::InvalidAttributeValue(tag, "id", szid);

					prec->SetItemList(pns);
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
				assert(tag.isend());
			}

			fem.AddDataRecord(prec);
		}
		else if (tag == "element_data")
		{
			sz = tag.AttributeValue("file", true);

			ElementDataRecord* prec = new ElementDataRecord(&fem, sz);
			const char* szdata = tag.AttributeValue("data");
			prec->Parse(szdata);

			const char* szname = tag.AttributeValue("name", true);
			if (szname != 0) prec->SetName(szname); else prec->SetName(szdata);

			sz = tag.AttributeValue("delim", true);
			if (sz != 0) prec->SetDelim(sz);

			sz = tag.AttributeValue("comments", true);
			if (sz != 0)
			{
				if      (strcmp(sz, "on") == 0) prec->SetComments(true);
				else if (strcmp(sz, "off") == 0) prec->SetComments(false); 
			}

			prec->SetItemList(tag.szvalue());

			fem.AddDataRecord(prec);
		}
		else if (tag == "rigid_body_data")
		{
			sz = tag.AttributeValue("file", true);

			ObjectDataRecord* prec = new ObjectDataRecord(&fem, sz);
			const char* szdata = tag.AttributeValue("data");
			prec->Parse(szdata);

			const char* szname = tag.AttributeValue("name", true);
			if (szname != 0) prec->SetName(szname); else prec->SetName(szdata);

			sz = tag.AttributeValue("delim", true);
			if (sz != 0) prec->SetDelim(sz);

			sz = tag.AttributeValue("comments", true);
			if (sz != 0)
			{
				if      (strcmp(sz, "on") == 0) prec->SetComments(true);
				else if (strcmp(sz, "off") == 0) prec->SetComments(false); 
			}

			prec->SetItemList(tag.szvalue());

			fem.AddDataRecord(prec);
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioOutputSection::ParsePlotfile(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	PlotFile* pplt = 0;
	const char* sz = tag.AttributeValue("type", true);
	if (sz)
	{
		if (strcmp(sz, "febio") == 0) pplt = new FEBioPlotFile(fem);
		else if (strcmp(sz, "lsdyna") == 0) pplt = new LSDYNAPlotFile;
		else throw XMLReader::InvalidAttributeValue(tag, "type", sz);
	}
	else pplt = new FEBioPlotFile(fem);
	fem.SetPlotFile(pplt);

	const char* szplt = tag.AttributeValue("file", true);
	if (szplt) m_pim->SetPlotfileName(szplt);

	if (dynamic_cast<LSDYNAPlotFile*>(pplt) && !tag.isleaf())
	{
		LSDYNAPlotFile& plt = *dynamic_cast<LSDYNAPlotFile*>(pplt);

		++tag;
		do
		{
			if (tag == "shell_strain") tag.value(plt.m_bsstrn);
			else if (tag == "map")
			{
				const char* szfield = tag.AttributeValue("field");
				const char* szval = tag.szvalue();
				if (strcmp(szfield, "displacement") == 0)
				{
					if (strcmp(szval, "DISPLACEMENT") == 0) plt.m_nfield[0] = PLOT_DISPLACEMENT;
					else throw XMLReader::InvalidValue(tag);
				}
				else if (strcmp(szfield, "velocity") == 0)
				{
					if (strcmp(szval, "NONE") == 0) plt.m_nfield[1] = PLOT_NONE;
					else if (strcmp(szval, "VELOCITY") == 0) plt.m_nfield[1] = PLOT_VELOCITY;
					else if (strcmp(szval, "FLUID_FLUX") == 0) plt.m_nfield[1] = PLOT_FLUID_FLUX;
					else if (strcmp(szval, "CONTACT_TRACTION") == 0) plt.m_nfield[1] = PLOT_CONTACT_TRACTION;
					else if (strcmp(szval, "REACTION_FORCE") == 0) plt.m_nfield[1] = PLOT_REACTION_FORCE;
					else if (strcmp(szval, "MATERIAL_FIBER") == 0) plt.m_nfield[1] = PLOT_MATERIAL_FIBER;
					else throw XMLReader::InvalidValue(tag);
				}
				else if (strcmp(szfield, "acceleration") == 0)
				{
					if (strcmp(szval, "NONE") == 0) plt.m_nfield[2] = PLOT_NONE;
					else if (strcmp(szval, "ACCELERATION") == 0) plt.m_nfield[2] = PLOT_ACCELERATION;
					else if (strcmp(szval, "FLUID_FLUX") == 0) plt.m_nfield[2] = PLOT_FLUID_FLUX;
					else if (strcmp(szval, "CONTACT_TRACTION") == 0) plt.m_nfield[2] = PLOT_CONTACT_TRACTION;
					else if (strcmp(szval, "REACTION_FORCE") == 0) plt.m_nfield[2] = PLOT_REACTION_FORCE;
					else if (strcmp(szval, "MATERIAL_FIBER") == 0) plt.m_nfield[2] = PLOT_MATERIAL_FIBER;
					else throw XMLReader::InvalidValue(tag);
				}
				else if (strcmp(szfield, "temperature") == 0)
				{
					if (strcmp(szval, "NONE") == 0) plt.m_nfield[3] = PLOT_NONE;
					else if (strcmp(szval, "FLUID_PRESSURE") == 0) plt.m_nfield[3] = PLOT_FLUID_PRESSURE;
					else if (strcmp(szval, "CONTACT_PRESSURE") == 0) plt.m_nfield[3] = PLOT_CONTACT_PRESSURE;
					else if (strcmp(szval, "CONTACT_GAP") == 0) plt.m_nfield[3] = PLOT_CONTACT_GAP;
					else throw XMLReader::InvalidValue(tag);
				}
				else if (strcmp(szfield, "plastic strain") == 0)
				{
					if      (strcmp(szval, "PLASTIC_STRAIN"  ) == 0) plt.m_nfield[4] = PLOT_PLASTIC_STRAIN;
					else if (strcmp(szval, "FIBER_STRAIN"    ) == 0) plt.m_nfield[4] = PLOT_FIBER_STRAIN;
					else if (strcmp(szval, "DEV_FIBER_STRAIN") == 0) plt.m_nfield[4] = PLOT_DEV_FIBER_STRAIN;
					else throw XMLReader::InvalidValue(tag);
				}
				else throw XMLReader::InvalidAttributeValue(tag, "field", szfield);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}
	else if (dynamic_cast<FEBioPlotFile*>(pplt))
	{
		// change the extension of the plot file to .xplt
		fem.SetPlotFileNameExtension(".xplt");

		if (!tag.isleaf())
		{
			FEBioPlotFile& plt = *dynamic_cast<FEBioPlotFile*>(pplt);
			++tag;
			do
			{
				if (tag == "var")
				{
					// get the variable name
					const char* szt = tag.AttributeValue("type");

					// get the item list
					vector<int> item;
					if (tag.isempty() == false)
					{
						// TODO: currently, this is only supported for domain variables, where
						//       the list is a list of materials
						vector<int> lmat;
						tag.value(lmat);
						// convert the material list to a domain list
						mesh.DomainListFromMaterial(lmat, item);
					}

					if (plt.AddVariable(szt, item) == false) throw XMLReader::InvalidAttributeValue(tag, "type", szt);
				}
				++tag;
			}
			while (!tag.isend());
		}
	}
}
