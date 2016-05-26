#include "stdafx.h"
#include "FEBioLoadsSection.h"
#include "FEBioMech/FEPointBodyForce.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
#include "FECore/BC.h"
#include "FECore/FESurfaceLoad.h"
#include "FECore/FEEdgeLoad.h"

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
		int nversion = m_pim->Version();
		if (nversion < 0x0200)
		{
			if      (tag == "force"      ) ParseNodalLoad(tag);
			else if (tag == "body_force" ) ParseBodyForce(tag);
			else if (tag == "heat_source") ParseBodyLoad (tag);
			else ParseSurfaceLoad(tag);
		}
		else if (nversion == 0x0200)
		{
			if      (tag == "nodal_load"  ) ParseNodalLoad    (tag);
			else if (tag == "surface_load") ParseSurfaceLoad20(tag);
			else if (tag == "edge_load"   ) ParseEdgeLoad     (tag);
			else if (tag == "body_load"   ) ParseBodyLoad20   (tag);
			else throw XMLReader::InvalidTag(tag);
		}
		else if (nversion > 0x0200)
		{
			if      (tag == "nodal_load"  ) ParseNodalLoad25  (tag);
			else if (tag == "surface_load") ParseSurfaceLoad25(tag);
			else if (tag == "edge_load"   ) ParseEdgeLoad25   (tag);
			else if (tag == "body_load"   ) ParseBodyLoad20   (tag);
			else throw XMLReader::InvalidTag(tag);
		}
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
				m_pim->value(tag, pf->m_a);
			}
			else if (tag == "node")
			{
				m_pim->value(tag, pf->m_inode); 
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
		FEBodyLoad* pf = fecore_new<FEBodyLoad>(FEBODYLOAD_ID, szt, &fem);
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
void FEBioLoadsSection::ParseBodyLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEBodyLoad* pbl = fecore_new<FEBodyLoad>(FEBODYLOAD_ID, tag.Name(), &fem);
	if (pbl == 0) throw XMLReader::InvalidTag(tag);
	FEParameterList& PL = pbl->GetParameterList();
	++tag;
	do
	{
		if (m_pim->ReadParameter(tag, PL) == false) throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
	fem.AddBodyLoad(pbl);
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseBodyLoad20(XMLTag& tag)
{
	const char* sztype = tag.AttributeValue("type");
	FEModel& fem = *GetFEModel();
	FEBodyLoad* pbl = fecore_new<FEBodyLoad>(FEBODYLOAD_ID, sztype, &fem);
	if (pbl == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
	FEParameterList& PL = pbl->GetParameterList();
	m_pim->ReadParameterList(tag, PL);
	fem.AddBodyLoad(pbl);
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseNodalLoad(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	DOFS& dofs = fem.GetDOFS();

	int nversion = m_pim->Version();
	if (nversion >= 0x0200)
	{
		// count how many nodal forces there are
		int ncnf = tag.children();

		// get the bc
		const char* sz = tag.AttributeValue("bc");
		int bc = dofs.GetDOF(sz);
		if (bc == -1) throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

		// get the load curve
		sz = tag.AttributeValue("lc");
		int lc = atoi(sz)-1;

		// see if there is a set defined
		const char* szset = tag.AttributeValue("set", true);
		if (szset)
		{
			// make sure this is a leaf tag
			if (tag.isleaf() == false) throw XMLReader::InvalidValue(tag);

			// find the node set
			FENodeSet* pns = mesh.FindNodeSet(szset);
			if (pns == 0) throw XMLReader::InvalidAttributeValue(tag, "set", szset);

			// see if the scale attribute is defined
			double scale = 1.0;
			tag.AttributeValue("scale", scale, true);

			FENodeSet& ns = *pns;
			int N = ns.size();
			for (int i=0; i<N; ++i)
			{
				int n = ns[i];
				// create new nodal force
				FENodalLoad* pfc = dynamic_cast<FENodalLoad*>(fecore_new<FEBoundaryCondition>(FEBC_ID, "nodal load", &fem));
				pfc->SetDOF(bc);
				pfc->SetLoad(scale, lc);
				pfc->AddNode(n);
				fem.AddNodalLoad(pfc);

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					GetStep()->AddModelComponent(pfc);
					pfc->Deactivate();
				}
			}
		}
		else
		{
			// read the prescribed data
			++tag;
			for (int i=0; i<ncnf; ++i)
			{
				// get the nodal ID
				int n = atoi(tag.AttributeValue("id"))-1;

				// get the load scale factor
				double scale = 0.0;
				m_pim->value(tag, scale);

				// create new nodal force
				FENodalLoad* pfc = dynamic_cast<FENodalLoad*>(fecore_new<FEBoundaryCondition>(FEBC_ID, "nodal load", &fem));
				pfc->SetDOF(bc);
				pfc->SetLoad(scale, lc);
				pfc->AddNode(n);
				fem.AddNodalLoad(pfc);

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					GetStep()->AddModelComponent(pfc);
					pfc->Deactivate();
				}

				++tag;
			}
		}
	}
	else
	{
		// count how many nodal forces there are
		int ncnf = tag.children();

		// read the prescribed data
		++tag;
		for (int i=0; i<ncnf; ++i)
		{
			int n = atoi(tag.AttributeValue("id"))-1, bc;
			const char* sz = tag.AttributeValue("bc");

			bc = dofs.GetDOF(sz);
			if (bc == -1) throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

			sz = tag.AttributeValue("lc");
			int lc = atoi(sz) - 1;

			double scale = 0.0;
			m_pim->value(tag, scale);

			FENodalLoad* pfc = dynamic_cast<FENodalLoad*>(fecore_new<FEBoundaryCondition>(FEBC_ID, "nodal load", &fem));
			pfc->SetDOF(bc);
			pfc->SetLoad(scale, lc);
			pfc->AddNode(n);
			fem.AddNodalLoad(pfc);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddModelComponent(pfc);
				pfc->Deactivate();
			}

			++tag;
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseNodalLoad25(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	DOFS& dofs = fem.GetDOFS();

	// get the bc
	const char* sz = tag.AttributeValue("bc");
	int bc = dofs.GetDOF(sz);
	if (bc == -1) throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

	// get the node set
	const char* szset = tag.AttributeValue("node_set");
	FENodeSet* nodeSet = mesh.FindNodeSet(szset);
	if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

	// create nodal load
	FENodalLoad* pfc = dynamic_cast<FENodalLoad*>(fecore_new<FEBoundaryCondition>(FEBC_ID, "nodal load", &fem));
	pfc->SetDOF(bc);
	pfc->AddNodes(*nodeSet);

	// add nodal load to model
	fem.AddNodalLoad(pfc);
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddModelComponent(pfc);
		pfc->Deactivate();
	}

	// read parameters
	FEParameterList& pl = pfc->GetParameterList();
	m_pim->ReadParameterList(tag, pl);
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseSurfaceLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// count how many pressure cards there are
	int npr = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(npr);
	fem.GetMesh().AddSurface(psurf);

	// create surface load
	FESurfaceLoad* ps = fecore_new<FESurfaceLoad>(FESURFACELOAD_ID, tag.Name(), &fem);
	if (ps == 0) throw XMLReader::InvalidTag(tag);

	// parse attributes
	for (int i=0; i<tag.m_natt; ++i)
	{
		XMLAtt& att = tag.m_att[i];
		if (ps->SetAttribute(att.m_szatt, att.m_szatv) == false) throw XMLReader::InvalidAttributeValue(tag, att.m_szatt, att.m_szatv);
	}

	// read the pressure data
	++tag;
	int nf[FEElement::MAX_NODES ], N;
	for (int i=0; i<npr; ++i)
	{
		FESurfaceElement& el = psurf->Element(i);

		for (int j=0; j<tag.m_natt; ++j)
		{
			XMLAtt& att = tag.m_att[j];
			if (ps->SetFacetAttribute(i, att.m_szatt, att.m_szatv) == false) throw XMLReader::InvalidAttributeValue(tag, att.m_szatt, att.m_szatv);
		}

		if      (tag == "quad4") el.SetType(FE_QUAD4G4);
		else if (tag == "tri3" ) el.SetType(m_pim->m_ntri3);
		else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
		else if (tag == "tri7" ) el.SetType(m_pim->m_ntri7);
		else if (tag == "quad8") el.SetType(FE_QUAD8G9);
		else if (tag == "quad9") el.SetType(FE_QUAD9G9);
		else throw XMLReader::InvalidTag(tag);

		N = el.Nodes();
		tag.value(nf, N);
		for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;

		++tag;
	}

	ps->SetSurface(psurf);

	// add surface load to model
	fem.AddSurfaceLoad(ps);

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddModelComponent(ps);
		ps->Deactivate();
	}
}


//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseSurfaceLoad20(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// create surface load
	const char* sztype = tag.AttributeValue("type");
	FESurfaceLoad* psl = fecore_new<FESurfaceLoad>(FESURFACELOAD_ID, sztype, &fem);
	if (psl == 0) throw XMLReader::InvalidTag(tag);

	// read name attribute
	const char* szname = tag.AttributeValue("name", true);
	if (szname) psl->SetName(szname);

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	fem.GetMesh().AddSurface(psurf);

	// read the parameters
	FEParameterList& pl = psl->GetParameterList();

	// read the pressure data
	++tag;
	do
	{
		if (m_pim->ReadParameter(tag, pl) == false)
		{
			if (tag == "surface")
			{
				// see if the surface is referenced by a set of defined explicitly
				const char* szset = tag.AttributeValue("set", true);
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
						if (m_pim->BuildSurface(*psurf, *ps) == false) throw XMLReader::InvalidTag(tag);
						psl->SetSurface(psurf);
					}
					else throw XMLReader::InvalidAttributeValue(tag, "set", szset);
				}
				else
				{
					// count how many pressure cards there are
					int npr = tag.children();
					psurf->create(npr);
					psl->SetSurface(psurf);

					++tag;
					int nf[FEElement::MAX_NODES ], N;
					for (int i=0; i<npr; ++i)
					{
						FESurfaceElement& el = psurf->Element(i);

						for (int j=0; j<tag.m_natt; ++j)
						{
							XMLAtt& att = tag.m_att[j];
							if (psl->SetFacetAttribute(i, att.m_szatt, att.m_szatv) == false) throw XMLReader::InvalidAttributeValue(tag, att.m_szatt, att.m_szatv);
						}

						if      (tag == "quad4") el.SetType(FE_QUAD4G4);
						else if (tag == "tri3" ) el.SetType(m_pim->m_ntri3);
						else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
						else if (tag == "tri7" ) el.SetType(m_pim->m_ntri7);
						else if (tag == "quad8") el.SetType(FE_QUAD8G9);
						else if (tag == "quad9") el.SetType(FE_QUAD9G9);
						else throw XMLReader::InvalidTag(tag);

						N = el.Nodes();
						tag.value(nf, N);
						for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;

						++tag;
					}
				}
			}
			else throw XMLReader::InvalidTag(tag);
		}
		++tag;
	}
	while (!tag.isend());

	// add surface load to model
	fem.AddSurfaceLoad(psl);

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddModelComponent(psl);
		psl->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseSurfaceLoad25(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// create surface load
	const char* sztype = tag.AttributeValue("type");
	FESurfaceLoad* psl = fecore_new<FESurfaceLoad>(FESURFACELOAD_ID, sztype, &fem);
	if (psl == 0) throw XMLReader::InvalidTag(tag);

	// read name attribute
	const char* szname = tag.AttributeValue("name", true);
	if (szname) psl->SetName(szname);

	// get the surface
	const char* szset = tag.AttributeValue("surface");
	FEFacetSet* pface = mesh.FindFacetSet(szset);
	if (pface == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", szset);

	FESurface* psurf = new FESurface(&mesh);
	m_pim->BuildSurface(*psurf, *pface);
	
	mesh.AddSurface(psurf);
	psl->SetSurface(psurf);

	// read the parameters
	FEParameterList& pl = psl->GetParameterList();
	m_pim->ReadParameterList(tag, pl);

	// add surface load to model
	fem.AddSurfaceLoad(psl);

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddModelComponent(psl);
		psl->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseEdgeLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// create edge load
	const char* sztype = tag.AttributeValue("type");
	FEEdgeLoad* pel = fecore_new<FEEdgeLoad>(FEEDGELOAD_ID, sztype, &fem);
	if (pel == 0) throw XMLReader::InvalidTag(tag);

	// create a new edge
	FEEdge* pedge = new FEEdge(&fem.GetMesh());
	fem.GetMesh().AddEdge(pedge);
	pel->SetEdge(pedge);

	// read the parameters
	FEParameterList& pl = pel->GetParameterList();

	// read the pressure data
	++tag;
	do
	{
		if (m_pim->ReadParameter(tag, pl) == false)
		{
			if (tag == "edge")
			{
				// see if the surface is referenced by a set of defined explicitly
				const char* szset = tag.AttributeValue("set", true);
				if (szset)
				{
					// make sure this tag does not have any children
					if (!tag.isleaf()) throw XMLReader::InvalidTag(tag);

					// see if we can find the facet set
					FEMesh& m = GetFEModel()->GetMesh();
					FESegmentSet* ps = m.FindSegmentSet(szset);

					// create an edge from the segment set
					if (ps)
					{
						if (BuildEdge(*pedge, *ps) == false) throw XMLReader::InvalidTag(tag);
						pel->Create(pedge->Elements());
					}
					else throw XMLReader::InvalidAttributeValue(tag, "set", szset);
				}
				else
				{
					// count how many load cards there are
					int npr = tag.children();
					pedge->create(npr);
					pel->Create(npr);

					++tag;
					int nf[FEElement::MAX_NODES ], N;
					for (int i=0; i<npr; ++i)
					{
						FELineElement& el = pedge->Element(i);

						for (int j=0; j<tag.m_natt; ++j)
						{
							XMLAtt& att = tag.m_att[j];
							if (pel->SetFacetAttribute(i, att.m_szatt, att.m_szatv) == false) throw XMLReader::InvalidAttributeValue(tag, att.m_szatt, att.m_szatv);
						}

						if      (tag == "line2") el.SetType(FE_LINE2G1);
						else throw XMLReader::InvalidTag(tag);

						N = el.Nodes();
						tag.value(nf, N);
						for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;

						++tag;
					}
				}
			}
			else throw XMLReader::InvalidTag(tag);
		}
		++tag;
	}
	while (!tag.isend());

	// add edge load to model
	fem.AddEdgeLoad(pel);

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddModelComponent(pel);
		pel->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseEdgeLoad25(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// create edge load
	const char* sztype = tag.AttributeValue("type");
	FEEdgeLoad* pel = fecore_new<FEEdgeLoad>(FEEDGELOAD_ID, sztype, &fem);
	if (pel == 0) throw XMLReader::InvalidTag(tag);

	// create a new edge
	FEEdge* pedge = new FEEdge(&fem.GetMesh());
	mesh.AddEdge(pedge);
	pel->SetEdge(pedge);

	// get the segment set
	const char* szedge = tag.AttributeValue("edge");
	FESegmentSet* pset = mesh.FindSegmentSet(szedge);
	if (pset == 0) throw XMLReader::InvalidAttributeValue(tag, "edge", szedge);
	if (BuildEdge(*pedge, *pset) == false) throw XMLReader::InvalidTag(tag);

	// read the parameters
	FEParameterList& pl = pel->GetParameterList();
	m_pim->ReadParameterList(tag, pl);

	// add edge load to model
	fem.AddEdgeLoad(pel);

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddModelComponent(pel);
		pel->Deactivate();
	}
}

//-----------------------------------------------------------------------------
bool FEBioLoadsSection::BuildEdge(FEEdge& e, FESegmentSet& es)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();
	int NN = m.Nodes();

	// count nr of segments
	int nsegs = es.Segments();

	// allocate storage for faces
	e.create(nsegs);

	// read segments
	for (int i=0; i<nsegs; ++i)
	{
		FELineElement& el = e.Element(i);
		FESegmentSet::SEGMENT& si = es.Segment(i);

		if      (si.ntype == 2) el.SetType(FE_LINE2G1);
		else return false;

		int N = el.Nodes(); assert(N == si.ntype);
		for (int j=0; j<N; ++j) el.m_node[j] = si.node[j];
	}

	// copy the name
	e.SetName(es.GetName());

	return true;
}
