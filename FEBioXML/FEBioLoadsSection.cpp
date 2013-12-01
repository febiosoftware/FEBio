#include "stdafx.h"
#include "FEBioLoadsSection.h"
#include "FEBioMech/FEPointBodyForce.h"
#include "FEBioHeat/FEHeatSource.h"
#include "FECore/FEModel.h"

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
		if      (tag == "force"              ) ParseBCForce    (tag);
		else if (tag == "body_force"         ) ParseBodyForce  (tag);
		else if (tag == "heat_source"        ) ParseHeatSource (tag);
		else ParseSurfaceLoad(tag);
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


//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseBCForce(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();

	int nversion = m_pim->Version();
	if (nversion >= 0x0200)
	{
		// count how many nodal forces there are
		int ncnf = tag.children();

		// get the bc
		int bc = -1;
		const char* sz = tag.AttributeValue("bc");

		if      (strcmp(sz, "x") == 0) bc = DOF_X;
		else if (strcmp(sz, "y") == 0) bc = DOF_Y;
		else if (strcmp(sz, "z") == 0) bc = DOF_Z;
		else if (strcmp(sz, "p") == 0) bc = DOF_P;
		else if (strcmp(sz, "t") == 0) bc = DOF_T;
		else if (strcmp(sz, "c") == 0) bc = DOF_C;
		else if (strcmp(sz, "c1") == 0) bc = DOF_C;
		else if (strcmp(sz, "c2") == 0) bc = DOF_C+1;
		else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

		// read the prescribed data
		++tag;
		for (int i=0; i<ncnf; ++i)
		{
			// get the nodal ID
			int n = atoi(tag.AttributeValue("id"))-1;

			// get the load curve
			sz = tag.AttributeValue("lc");
			int lc = atoi(sz)-1;

			// create new nodal force
			FENodalForce* pfc = new FENodalForce;
			pfc->node = n;
			pfc->bc = bc;
			pfc->lc = lc;
			tag.value(pfc->s);
			fem.AddNodalLoad(pfc);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddBoundaryCondition(pfc);
				pfc->Deactivate();
			}

			++tag;
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

			if      (strcmp(sz, "x") == 0) bc = 0;
			else if (strcmp(sz, "y") == 0) bc = 1;
			else if (strcmp(sz, "z") == 0) bc = 2;
			else if (strcmp(sz, "p") == 0) bc = 6;
			else if (strcmp(sz, "t") == 0) bc = 10;
			else if (strcmp(sz, "c") == 0) bc = 11;
			else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

			sz = tag.AttributeValue("lc");
			int lc = atoi(sz) - 1;

			FENodalForce* pfc = new FENodalForce;
			pfc->node = n;
			pfc->bc = bc;
			pfc->lc = lc;
			tag.value(pfc->s);
			fem.AddNodalLoad(pfc);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddBoundaryCondition(pfc);
				pfc->Deactivate();
			}

			++tag;
		}
	}
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
	FEBioKernel& febio = FEBioKernel::GetInstance();
	FESurfaceLoad* ps = febio.Create<FESurfaceLoad>(tag.Name(), &fem);
	if (ps == 0) throw XMLReader::InvalidTag(tag);

	ps->Create(npr);
	ps->SetSurface(psurf);

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
		else if (tag == "quad8") el.SetType(FE_QUAD8G9);
		else throw XMLReader::InvalidTag(tag);

		N = el.Nodes();
		tag.value(nf, N);
		for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;

		++tag;
	}

	// add surface load to model
	fem.AddSurfaceLoad(ps);

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(ps);
		ps->Deactivate();
	}
}
