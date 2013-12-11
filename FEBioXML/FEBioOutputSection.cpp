#include "stdafx.h"
#include "FEBioOutputSection.h"
#include "FECore/NodeDataRecord.h"
#include "FECore/ElementDataRecord.h"
#include "FECore/ObjectDataRecord.h"
#include "FECore/FEModel.h"

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

	// get the plot file type. Must be "febio"!
	const char* sz = tag.AttributeValue("type", true);
	if (sz)
	{
		if (strcmp(sz, "febio") != 0) throw XMLReader::InvalidAttributeValue(tag, "type", sz);
	}
	else sz = "febio";
	strcpy(m_pim->m_szplot_type, sz);

	// get the optional plot file name
	const char* szplt = tag.AttributeValue("file", true);
	if (szplt) m_pim->SetPlotfileName(szplt);

	// read and store the plot variables
	if (!tag.isleaf())
	{
		++tag;
		do
		{
			if (tag == "var")
			{
				// get the variable name
				const char* szt = tag.AttributeValue("type");

				// get the item list
				vector<int> item;
				if (tag.isempty() == false) tag.value(item);

				// Add the plot variable
				m_pim->AddPlotVariable(szt, item);
			}
			++tag;
		}
		while (!tag.isend());
	}
}
