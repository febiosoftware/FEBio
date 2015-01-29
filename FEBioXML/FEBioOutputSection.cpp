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

	// see if the log file has any attributes
	const char* szlog = tag.AttributeValue("file", true);
	if (szlog) m_pim->SetLogfileName(szlog);

	if (tag.isleaf()) return;
	++tag;
	do
	{
		if (tag == "node_data")
		{
			const char* sz = tag.AttributeValue("file", true);

			// if we have a path, prepend the path's name
			const char* szpath = m_pim->m_szpath;
			char szfile[1024] = {0};
			if (szpath && szpath[0])
			{
				sprintf(szfile, "%s%s", szpath, sz);
			}
			else strcpy(szfile, sz);

			NodeDataRecord* prec = new NodeDataRecord(&fem, szfile);
			const char* szdata = tag.AttributeValue("data");
			prec->Parse(szdata);

			const char* szname = tag.AttributeValue("name", true);
			if (szname != 0) prec->SetName(szname); else prec->SetName(szdata);

			sz = tag.AttributeValue("delim", true);
			if (sz != 0) prec->SetDelim(sz);

			sz = tag.AttributeValue("format", true);
			if (sz!=0) prec->SetFormat(sz);

			sz = tag.AttributeValue("comments", true);
			if (sz != 0)
			{
				if      (strcmp(sz, "on") == 0) prec->SetComments(true);
				else if (strcmp(sz, "off") == 0) prec->SetComments(false); 
			}

			sz = tag.AttributeValue("set", true);
			if (sz)
			{
				FENodeSet* pns = mesh.FindNodeSet(sz);
				if (pns == 0) throw XMLReader::InvalidAttributeValue(tag, "set", sz);
				prec->SetItemList(pns);
			}
			else prec->DataRecord::SetItemList(tag.szvalue());

			fem.AddDataRecord(prec);
		}
		else if (tag == "element_data")
		{
			const char* sz = tag.AttributeValue("file", true);

			// if we have a path, prepend the path's name
			const char* szpath = m_pim->m_szpath;
			char szfile[1024] = {0};
			if (szpath && szpath[0])
			{
				sprintf(szfile, "%s%s", szpath, sz);
			}
			else strcpy(szfile, sz);

			ElementDataRecord* prec = new ElementDataRecord(&fem, szfile);
			const char* szdata = tag.AttributeValue("data");
			prec->Parse(szdata);

			const char* szname = tag.AttributeValue("name", true);
			if (szname != 0) prec->SetName(szname); else prec->SetName(szdata);

			sz = tag.AttributeValue("delim", true);
			if (sz != 0) prec->SetDelim(sz);

			sz = tag.AttributeValue("format", true);
			if (sz!=0) prec->SetFormat(sz);

			sz = tag.AttributeValue("comments", true);
			if (sz != 0)
			{
				if      (strcmp(sz, "on") == 0) prec->SetComments(true);
				else if (strcmp(sz, "off") == 0) prec->SetComments(false); 
			}

			sz = tag.AttributeValue("elset", true);
			if (sz)
			{
				FEElementSet* pes = mesh.FindElementSet(sz);
				if (pes == 0) throw XMLReader::InvalidAttributeValue(tag, "elset", sz);
				prec->SetItemList(pes);
			}
			else prec->DataRecord::SetItemList(tag.szvalue());

			fem.AddDataRecord(prec);
		}
		else if (tag == "rigid_body_data")
		{
			const char* sz = tag.AttributeValue("file", true);

			// if we have a path, prepend the path's name
			const char* szpath = m_pim->m_szpath;
			char szfile[1024] = {0};
			if (szpath && szpath[0])
			{
				sprintf(szfile, "%s%s", szpath, sz);
			}
			else strcpy(szfile, sz);

			ObjectDataRecord* prec = new ObjectDataRecord(&fem, szfile);
			const char* szdata = tag.AttributeValue("data");
			prec->Parse(szdata);

			const char* szname = tag.AttributeValue("name", true);
			if (szname != 0) prec->SetName(szname); else prec->SetName(szdata);

			sz = tag.AttributeValue("delim", true);
			if (sz != 0) prec->SetDelim(sz);

			sz = tag.AttributeValue("format", true);
			if (sz!=0) prec->SetFormat(sz);

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
			else if (tag=="compression")
			{
				int ncomp;
				tag.value(ncomp);
				m_pim->SetPlotCompression(ncomp);
			}
			++tag;
		}
		while (!tag.isend());
	}
}
