#include "stdafx.h"
#include "FEBioOutputSection.h"
#include "FECore/NodeDataRecord.h"
#include "FECore/ElementDataRecord.h"
#include "FECore/ObjectDataRecord.h"
#include "FECore/NLConstraintDataRecord.h"
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

			NodeDataRecord* prec = 0;
			if (sz)
			{
				// if we have a path, prepend the path's name
				const char* szpath = m_pim->m_szpath;
				char szfile[1024] = {0};
				if (szpath && szpath[0])
				{
					sprintf(szfile, "%s%s", szpath, sz);
				}
				else strcpy(szfile, sz);

				prec = new NodeDataRecord(&fem, szfile);
			}
			else prec = new NodeDataRecord(&fem, 0);

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

			m_pim->AddDataRecord(prec);
		}
		else if (tag == "element_data")
		{
			const char* sz = tag.AttributeValue("file", true);

			ElementDataRecord* prec = 0;
			if (sz)
			{
				// if we have a path, prepend the path's name
				const char* szpath = m_pim->m_szpath;
				char szfile[1024] = {0};
				if (szpath && szpath[0])
				{
					sprintf(szfile, "%s%s", szpath, sz);
				}
				else strcpy(szfile, sz);
				prec = new ElementDataRecord(&fem, szfile);
			}
			else prec = new ElementDataRecord(&fem, 0);

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

			const char* sztmp = "elset";
			if (m_pim->Version() >= 0x0205) sztmp = "elem_set";

			sz = tag.AttributeValue(sztmp, true);
			if (sz)
			{
				FEElementSet* pes = mesh.FindElementSet(sz);
				if (pes == 0) throw XMLReader::InvalidAttributeValue(tag, sztmp, sz);
				prec->SetItemList(pes);
			}
			else prec->DataRecord::SetItemList(tag.szvalue());

			m_pim->AddDataRecord(prec);
		}
		else if (tag == "rigid_body_data")
		{
			const char* sz = tag.AttributeValue("file", true);

			ObjectDataRecord* prec = 0;
			if (sz)
			{
				// if we have a path, prepend the path's name
				const char* szpath = m_pim->m_szpath;
				char szfile[1024] = {0};
				if (szpath && szpath[0])
				{
					sprintf(szfile, "%s%s", szpath, sz);
				}
				else strcpy(szfile, sz);
				prec = new ObjectDataRecord(&fem, szfile);
			}
			else prec = new ObjectDataRecord(&fem, 0);

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

			m_pim->AddDataRecord(prec);
		}
        else if (tag == "rigid_connector_data")
        {
            const char* sz = tag.AttributeValue("file", true);
            
            NLConstraintDataRecord* prec = 0;
            if (sz)
            {
                // if we have a path, prepend the path's name
                const char* szpath = m_pim->m_szpath;
                char szfile[1024] = {0};
                if (szpath && szpath[0])
                {
                    sprintf(szfile, "%s%s", szpath, sz);
                }
                else strcpy(szfile, sz);
                prec = new NLConstraintDataRecord(&fem, szfile);
            }
            else prec = new NLConstraintDataRecord(&fem, 0);
            
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
            
            m_pim->AddDataRecord(prec);
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
		if ((strcmp(sz, "febio") != 0) && (strcmp(sz, "febio2") != 0)) throw XMLReader::InvalidAttributeValue(tag, "type", sz);
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

                // see if a surface is referenced
                const char* szset = tag.AttributeValue("surface", true);
                if (szset)
                {
                    // make sure this tag does not have any children
                    if (!tag.isleaf()) throw XMLReader::InvalidTag(tag);
                    
                    // see if we can find the facet set
                    FEMesh& m = GetFEModel()->GetMesh();
                    FEFacetSet* ps = 0;
                    for (int i=0; i<m.FacetSets(); ++i)
                    {
                        FEFacetSet& fi = m.FacetSet(i);
                        if (strcmp(fi.GetName(), szset) == 0)
                        {
                            ps = &fi;
                            break;
                        }
                    }
                    
                    // create a surface from the facet set
                    if (ps)
                    {
                        // create a new surface
                        FESurface* psurf = new FESurface(&fem.GetMesh());
                        fem.GetMesh().AddSurface(psurf);
                        if (BuildSurface(*psurf, *ps) == false) throw XMLReader::InvalidTag(tag);
                        // Add the plot variable
                        const char* szd = psurf->GetName();
                        m_pim->AddPlotVariable(szt, item, szd);
                    }
                    else throw XMLReader::InvalidAttributeValue(tag, "set", szset);
                }
                else
                {
                    // Add the plot variable
                    m_pim->AddPlotVariable(szt, item);
                }
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


//-----------------------------------------------------------------------------
bool FEBioOutputSection::BuildSurface(FESurface& s, FEFacetSet& fs)
{
    FEModel& fem = *GetFEModel();
    FEMesh& m = fem.GetMesh();
    int NN = m.Nodes();
    
    // count nr of faces
    int faces = fs.Faces();
    
    // allocate storage for faces
    s.Create(faces);
    
    // read faces
    for (int i=0; i<faces; ++i)
    {
        FESurfaceElement& el = s.Element(i);
        FEFacetSet::FACET& fi = fs.Face(i);
        
        if      (fi.ntype == 4) el.SetType(FE_QUAD4G4);
        else if (fi.ntype == 3) el.SetType(m_pim->m_ntri3);
        else if (fi.ntype == 6) el.SetType(m_pim->m_ntri6);
        else if (fi.ntype == 7) el.SetType(m_pim->m_ntri7);
        else if (fi.ntype == 8) el.SetType(FE_QUAD8G9);
        else if (fi.ntype == 9) el.SetType(FE_QUAD9G9);
        else return false;
        
        int N = el.Nodes(); assert(N == fi.ntype);
        for (int j=0; j<N; ++j) el.m_node[j] = fi.node[j];
    }
    
    // copy the name
    s.SetName(fs.GetName());
    
    return true;
}
