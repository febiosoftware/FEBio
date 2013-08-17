#include "stdafx.h"
#include "FEBioBoundarySection.h"
#include "FEBioMech/FEPressureLoad.h"
#include "FEBioMech/FETractionLoad.h"
#include "FEBioMix/FEPoroTraction.h"
#include "FEBioMix/FEFluidFlux.h"
#include "FEBioMix/FESoluteFlux.h"
#include "FEBioMech/FEDiscreteMaterial.h"
#include "FEBioMech/FEDiscreteSpringDomain.h"
#include "FEBioMech/FESlidingInterface.h"
#include "FEBioMix/FESlidingInterface2.h"
#include "FEBioMix/FESlidingInterface3.h"
#include "FEBioMix/FETiedBiphasicInterface.h"
#include "FEBioMech/FETiedInterface.h"
#include "FEBioMech/FEFacet2FacetSliding.h"
#include "FEBioMech/FEFacet2FacetTied.h"
#include "FEBioMech/FESlidingInterfaceBW.h"
#include "FEBioMech/FEPeriodicBoundary.h"
#include "FEBioMech/FESurfaceConstraint.h"
#include "FEBioMech/FERigidWallInterface.h"
#include "FEBioMech/FEAugLagLinearConstraint.h"
#include "FEBioMech/FERigidJoint.h"
#include "FEBioHeat/FEHeatFlux.h"
#include "FEBioHeat/FEConvectiveHeatFlux.h"

//-----------------------------------------------------------------------------
//!  Parses the boundary section from the xml file
//!
void FEBioBoundarySection::Parse(XMLTag& tag)
{
	// make sure this tag has children
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "fix"                  ) ParseBCFix               (tag);
		else if (tag == "prescribe"            ) ParseBCPrescribe         (tag);
		else if (tag == "contact"              ) ParseContactSection      (tag);
		else if (tag == "linear_constraint"    ) ParseConstraints         (tag);
		else if (tag == "spring"               ) ParseSpringSection       (tag);
		else if (m_pim->Version() < 0x0102)
		{
			// As of version 1.2 and up, the following functions are handled by the Loads section
			if      (tag == "force"                ) ParseBCForce             (tag);
			else if (tag == "pressure"             ) ParseBCPressure          (tag);
			else if (tag == "traction"             ) ParseBCTraction          (tag);
			else if (tag == "normal_traction"      ) ParseBCPoroNormalTraction(tag);
			else if (tag == "fluidflux"            ) ParseBCFluidFlux         (tag);
			else if (tag == "soluteflux"           ) ParseBCSoluteFlux        (tag);
			else if (tag == "heatflux"             ) ParseBCHeatFlux          (tag);
			else throw XMLReader::InvalidTag(tag);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//---------------------------------------------------------------------------------
// parse a surface section for contact definitions
//
bool FEBioBoundarySection::ParseSurfaceSection(XMLTag &tag, FESurface& s, int nfmt, bool bnodal)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();
	int NN = m.Nodes();

	int N, nf[8];

	// count nr of faces
	int faces = tag.children();

	// allocate storage for faces
	s.create(faces);

	// read faces
	++tag;
	for (int i=0; i<faces; ++i)
	{
		FESurfaceElement& el = s.Element(i);

		// set the element type/integration rule
		if (bnodal)
		{
			if      (tag == "quad4") el.SetType(FE_QUAD4NI);
			else if (tag == "tri3" ) el.SetType(FE_TRI3NI );
			else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
			else throw XMLReader::InvalidTag(tag);
		}
		else
		{
			if      (tag == "quad4") el.SetType(FE_QUAD4G4);
			else if (tag == "tri3" ) el.SetType(m_pim->m_ntri3);
			else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
			else if (tag == "quad8") el.SetType(FE_QUAD8G9);
			else throw XMLReader::InvalidTag(tag);
		}

		N = el.Nodes();

		if (nfmt == 0)
		{
			tag.value(nf, N);
			for (int j=0; j<N; ++j) 
			{
				int nid = nf[j]-1;
				if ((nid<0)||(nid>= NN)) throw XMLReader::InvalidValue(tag);
				el.m_node[j] = nid;
			}
		}
		else if (nfmt == 1)
		{
			tag.value(nf, 2);
			FEElement* pe = m.FindElementFromID(nf[0]);
			if (pe)
			{
				int ne[8];
				int nn = m.GetFace(*pe, nf[1]-1, ne);
				if (nn != N) throw XMLReader::InvalidValue(tag);
				for (int j=0; j<N; ++j) el.m_node[j] = ne[j];
				el.m_nelem = nf[0];
			}
			else throw XMLReader::InvalidValue(tag);
		}

		++tag;
	}
	return true;
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseBCFix(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();

	// make sure this section does not appear in a step section
	if (m_pim->m_nsteps != 0) throw XMLReader::InvalidTag(tag);

	// see if a set is defined
	const char* szset = tag.AttributeValue("set", true);
	if (szset)
	{
		// read the set
		FEMesh& mesh = fem.GetMesh();
		FENodeSet* ps = mesh.FindNodeSet(szset);
		if (ps == 0) throw XMLReader::InvalidAttributeValue(tag, "set", szset);

		// get the bc attribute
		const char* sz = tag.AttributeValue("bc");

		// Make sure this is a leaf
		if (tag.isleaf() == false) throw XMLReader::InvalidValue(tag);

		// loop over all nodes in the nodeset
		FENodeSet& s = *ps;
		int N = s.size();
		for (int i=0; i<N; ++i)
		{
			FENode& node = mesh.Node(s[i]);
			if      (strcmp(sz, "x"  ) == 0) { node.m_ID[DOF_X] = -1; }
			else if (strcmp(sz, "y"  ) == 0) { node.m_ID[DOF_Y] = -1; }
			else if (strcmp(sz, "z"  ) == 0) { node.m_ID[DOF_Z] = -1; }
			else if (strcmp(sz, "xy" ) == 0) { node.m_ID[DOF_X] = node.m_ID[DOF_Y] = -1; }
			else if (strcmp(sz, "yz" ) == 0) { node.m_ID[DOF_Y] = node.m_ID[DOF_Z] = -1; }
			else if (strcmp(sz, "xz" ) == 0) { node.m_ID[DOF_X] = node.m_ID[DOF_Z] = -1; }
			else if (strcmp(sz, "xyz") == 0) { node.m_ID[DOF_X] = node.m_ID[DOF_Y] = node.m_ID[DOF_Z] = -1; }
			else if (strcmp(sz, "p"  ) == 0) { node.m_ID[DOF_P] = -1; }
			else if (strcmp(sz, "u"  ) == 0) { node.m_ID[DOF_U] = -1; }
			else if (strcmp(sz, "v"  ) == 0) { node.m_ID[DOF_V] = -1; }
			else if (strcmp(sz, "w"  ) == 0) { node.m_ID[DOF_W] = -1; }
			else if (strcmp(sz, "uv" ) == 0) { node.m_ID[DOF_U] = node.m_ID[DOF_V] = -1; }
			else if (strcmp(sz, "vw" ) == 0) { node.m_ID[DOF_V] = node.m_ID[DOF_W] = -1; }
			else if (strcmp(sz, "uw" ) == 0) { node.m_ID[DOF_U] = node.m_ID[DOF_W] = -1; }
			else if (strcmp(sz, "uvw") == 0) { node.m_ID[DOF_U] = node.m_ID[DOF_V] = node.m_ID[DOF_W] = -1; }
			else if (strcmp(sz, "t"  ) == 0) { node.m_ID[DOF_T] = -1; }
			else if (strcmp(sz, "c"  ) == 0) { node.m_ID[DOF_C] = -1; }
			else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);
		}
	}
	else
	{
		// Read the fixed nodes
		++tag;
		do
		{
			int n = atoi(tag.AttributeValue("id"))-1;
			FENode& node = fem.GetMesh().Node(n);
			const char* sz = tag.AttributeValue("bc");
			if      (strcmp(sz, "x") == 0) node.m_ID[DOF_X] = -1;
			else if (strcmp(sz, "y") == 0) node.m_ID[DOF_Y] = -1;
			else if (strcmp(sz, "z") == 0) node.m_ID[DOF_Z] = -1;
			else if (strcmp(sz, "xy") == 0) { node.m_ID[DOF_X] = node.m_ID[DOF_Y] = -1; }
			else if (strcmp(sz, "yz") == 0) { node.m_ID[DOF_Y] = node.m_ID[DOF_Z] = -1; }
			else if (strcmp(sz, "xz") == 0) { node.m_ID[DOF_X] = node.m_ID[DOF_Z] = -1; }
			else if (strcmp(sz, "xyz") == 0) { node.m_ID[DOF_X] = node.m_ID[DOF_Y] = node.m_ID[DOF_Z] = -1; }
			else if (strcmp(sz, "p") == 0) { node.m_ID[DOF_P] = -1; }
			else if (strcmp(sz, "u") == 0) node.m_ID[DOF_U] = -1;
			else if (strcmp(sz, "v") == 0) node.m_ID[DOF_V] = -1;
			else if (strcmp(sz, "w") == 0) node.m_ID[DOF_W] = -1;
			else if (strcmp(sz, "uv") == 0) { node.m_ID[DOF_U] = node.m_ID[DOF_V] = -1; }
			else if (strcmp(sz, "vw") == 0) { node.m_ID[DOF_V] = node.m_ID[DOF_W] = -1; }
			else if (strcmp(sz, "uw") == 0) { node.m_ID[DOF_U] = node.m_ID[DOF_W] = -1; }
			else if (strcmp(sz, "uvw") == 0) { node.m_ID[DOF_U] = node.m_ID[DOF_V] = node.m_ID[DOF_W] = -1; }
			else if (strcmp(sz, "t") == 0) node.m_ID[DOF_T] = -1;
			else if (strcmp(sz, "c") == 0) { node.m_ID[DOF_C] = -1; }
			else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);
			++tag;
		}
		while (!tag.isend());
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseBCPrescribe(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	int nversion = m_pim->Version();

	if (nversion >= 0x0200)
	{
		// count how many prescibed nodes there are
		int ndis = tag.children();

		// determine whether prescribed BC is relative or absolute
		bool br = false;
		const char* sztype = tag.AttributeValue("type",true);
		if (sztype && strcmp(sztype, "relative") == 0) br = true;

		// get the BC
		int bc = -1;
		const char* sz = tag.AttributeValue("bc");
		if      (strcmp(sz, "x") == 0) bc = DOF_X;
		else if (strcmp(sz, "y") == 0) bc = DOF_Y;
		else if (strcmp(sz, "z") == 0) bc = DOF_Z;
		else if (strcmp(sz, "u") == 0) bc = DOF_U;
		else if (strcmp(sz, "v") == 0) bc = DOF_V;
		else if (strcmp(sz, "w") == 0) bc = DOF_W;
		else if (strcmp(sz, "p") == 0) bc = DOF_P;
		else if (strcmp(sz, "t") == 0) bc = DOF_T; 
		else if (strcmp(sz, "c") == 0) bc = DOF_C;
		else if (strncmp(sz, "c", 1) == 0) bc = DOF_C + atoi(&sz[1]) - 1;
		else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

		// read the prescribed data
		++tag;
		for (int i=0; i<ndis; ++i)
		{
			// get the node ID
			int n = atoi(tag.AttributeValue("id"))-1;

			// get the load curve number
			sz = tag.AttributeValue("lc");
			int lc = atoi(sz) - 1;

			// create a new BC
			FEPrescribedBC* pdc = new FEPrescribedBC;
			pdc->node = n;
			pdc->bc = bc;
			pdc->lc = lc;
			tag.value(pdc->s);
			pdc->br = br;
			fem.AddPrescribedBC(pdc);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddBoundaryCondition(pdc);
				pdc->Deactivate();
			}
			++tag;
		}
	}
	else
	{
		// see if this tag defines a set
		const char* szset = tag.AttributeValue("set", true);
		if (szset)
		{
			// Find the set
			FENodeSet* ps = mesh.FindNodeSet(szset);
			if (ps == 0) throw XMLReader::InvalidAttributeValue(tag, "set", szset);

			// get the bc attribute
			const char* sz = tag.AttributeValue("bc");

			int bc;
			if      (strcmp(sz, "x") == 0) bc = DOF_X;
			else if (strcmp(sz, "y") == 0) bc = DOF_Y;
			else if (strcmp(sz, "z") == 0) bc = DOF_Z;
			else if (strcmp(sz, "u") == 0) bc = DOF_U;
			else if (strcmp(sz, "v") == 0) bc = DOF_V;
			else if (strcmp(sz, "w") == 0) bc = DOF_W;
			else if (strcmp(sz, "p") == 0) bc = DOF_P;
			else if (strcmp(sz, "t") == 0) bc = DOF_T; 
			else if (strcmp(sz, "c") == 0) bc = DOF_C;
			else if (strncmp(sz, "c", 1) == 0) bc = DOF_C + atoi(&sz[1]) - 1;
			else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

			// get the lc attribute
			sz = tag.AttributeValue("lc");
			int lc = atoi(sz);

			// make sure this tag is a leaf
			if (tag.isleaf() == false) throw XMLReader::InvalidValue(tag);

			// get the scale factor
			double s = 1;
			tag.value(s);

			// loop over all nodes in the nodeset
			FENodeSet& ns = *ps;
			int N = ns.size();
			for (int i=0; i<N; ++i)
			{
				FEPrescribedBC* pdc = new FEPrescribedBC;
				pdc->node = ns[i];
				pdc->bc = bc;
				pdc->lc = lc;
				pdc->s = s;
				fem.AddPrescribedBC(pdc);

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					GetStep()->AddBoundaryCondition(pdc);
					pdc->Deactivate();
				}
			}
		}
		else
		{
			// count how many prescibed nodes there are
			int ndis = tag.children();

			// determine whether prescribed BC is relative or absolute
			bool br = false;
			const char* sztype = tag.AttributeValue("type",true);
			if (sztype && strcmp(sztype, "relative") == 0) br = true;

			// read the prescribed data
			++tag;
			for (int i=0; i<ndis; ++i)
			{
				int n = atoi(tag.AttributeValue("id"))-1, bc;
				const char* sz = tag.AttributeValue("bc");

				if      (strcmp(sz, "x") == 0) bc = DOF_X;
				else if (strcmp(sz, "y") == 0) bc = DOF_Y;
				else if (strcmp(sz, "z") == 0) bc = DOF_Z;
				else if (strcmp(sz, "u") == 0) bc = DOF_U;
				else if (strcmp(sz, "v") == 0) bc = DOF_V;
				else if (strcmp(sz, "w") == 0) bc = DOF_W;
				else if (strcmp(sz, "p") == 0) bc = DOF_P;
				else if (strcmp(sz, "t") == 0) bc = DOF_T; 
				else if (strcmp(sz, "c") == 0) bc = DOF_C;
				else if (strcmp(sz, "c1") == 0) bc = DOF_C;
				else if (strncmp(sz, "c", 1) == 0) bc = DOF_C + atoi(&sz[1]) - 1;
				else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

				sz = tag.AttributeValue("lc");
				int lc = atoi(sz)-1;

				FEPrescribedBC* pdc = new FEPrescribedBC;
				pdc->node = n;
				pdc->bc = bc;
				pdc->lc = lc;
				tag.value(pdc->s);
				pdc->br = br;
				fem.AddPrescribedBC(pdc);

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					GetStep()->AddBoundaryCondition(pdc);
					pdc->Deactivate();
				}
				++tag;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseBCForce(XMLTag &tag)
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
void FEBioBoundarySection::ParseBCPressure(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	const char* sz;
	bool blinear = false;
	sz = tag.AttributeValue("type", true);
	if (sz)
	{
		if (strcmp(sz, "linear") == 0) blinear = true;
		else if (strcmp(sz, "nonlinear") == 0) blinear = false;
		else throw XMLReader::InvalidAttributeValue(tag, "type", sz);
	}

	// count how many pressure cards there are
	int npr = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(npr);
	fem.GetMesh().AddSurface(psurf);

	// allocate pressure data
	FEPressureLoad* ps = new FEPressureLoad(psurf, blinear);
	ps->create(npr);
	fem.AddSurfaceLoad(ps);

	// read the pressure data
	++tag;
	int nf[FEElement::MAX_NODES ], N;
	double s;
	for (int i=0; i<npr; ++i)
	{
		FEPressureLoad::LOAD& pc = ps->PressureLoad(i);
		FESurfaceElement& el = psurf->Element(i);

		sz = tag.AttributeValue("lc");
		pc.lc = atoi(sz)-1;

		s  = atof(tag.AttributeValue("scale"));
		pc.s[0] = pc.s[1] = pc.s[2] = pc.s[3] = s;
		pc.s[4] = pc.s[5] = pc.s[6] = pc.s[7] = s;

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

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(ps);
		ps->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseBCTraction(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();

	const char* sz;

	// count how many traction cards there are
	int ntc = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(ntc);
	fem.GetMesh().AddSurface(psurf);

	// allocate traction data
	FETractionLoad* pt = new FETractionLoad(psurf);
	fem.AddSurfaceLoad(pt);
	pt->create(ntc);

	// read the traction data
	++tag;
	int nf[FEElement::MAX_NODES], N;
	vec3d s;
	for (int i=0; i<ntc; ++i)
	{
		FETractionLoad::LOAD& tc = pt->TractionLoad(i);
		FESurfaceElement& el = psurf->Element(i);

		sz = tag.AttributeValue("lc");
		tc.lc = atoi(sz)-1;

		s.x  = atof(tag.AttributeValue("tx"));
		s.y  = atof(tag.AttributeValue("ty"));
		s.z  = atof(tag.AttributeValue("tz"));

		tc.s[0] = tc.s[1] = tc.s[2] = tc.s[3] = s;
		tc.s[4] = tc.s[5] = tc.s[6] = tc.s[7] = s;

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

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(pt);
		pt->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseBCPoroNormalTraction(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	
	const char* sz;
	bool blinear = false;
	sz = tag.AttributeValue("type", true);
	if (sz)
	{
		if (strcmp(sz, "linear") == 0) blinear = true;
		else if (strcmp(sz, "nonlinear") == 0) blinear = false;
		else throw XMLReader::InvalidAttributeValue(tag, "type", sz);
	}
	
	bool beffective = false;
	sz = tag.AttributeValue("traction", true);
	if (sz)
	{
		if (strcmp(sz, "effective") == 0) beffective = true;
		else if ((strcmp(sz, "total") == 0) || (strcmp(sz, "mixture") == 0)) beffective = false;
		else throw XMLReader::InvalidAttributeValue(tag, "traction", sz);
	}
	
	// count how many normal traction cards there are
	int npr = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(npr);
	fem.GetMesh().AddSurface(psurf);
	
	// allocate normal traction data
	FEPoroNormalTraction* ps = new FEPoroNormalTraction(psurf, blinear, beffective);
	ps->create(npr);
	fem.AddSurfaceLoad(ps);
	
	// read the normal traction data
	++tag;
	int nf[4], N;
	double s;
	for (int i=0; i<npr; ++i)
	{
		FEPoroNormalTraction::LOAD& pc = ps->NormalTraction(i);
		FESurfaceElement& el = psurf->Element(i);
		
		sz = tag.AttributeValue("lc");
		pc.lc = atoi(sz)-1;
		
		s  = atof(tag.AttributeValue("scale"));
		pc.s[0] = pc.s[1] = pc.s[2] = pc.s[3] = s;
		pc.s[4] = pc.s[5] = pc.s[6] = pc.s[7] = s;
		
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

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(ps);
		ps->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseBCFluidFlux(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();

	const char* sz;
	bool blinear = false;
	sz = tag.AttributeValue("type", true);
	if (sz)
	{
		if (strcmp(sz, "linear") == 0) blinear = true;
		else if (strcmp(sz, "nonlinear") == 0) blinear = false;
		else throw XMLReader::InvalidAttributeValue(tag, "type", sz);
	}
	
	bool bmixture = false;
	sz = tag.AttributeValue("flux", true);
	if (sz)
	{
		if (strcmp(sz, "mixture") == 0) bmixture = true;
		else if (strcmp(sz, "fluid") == 0) bmixture = false;
		else throw XMLReader::InvalidAttributeValue(tag, "flux", sz);
	}
	
	// count how many fluid flux cards there are
	int nfr = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(nfr);
	fem.GetMesh().AddSurface(psurf);
	
	// allocate fluid flux data
	FEFluidFlux* pfs = new FEFluidFlux(psurf, blinear, bmixture);
	pfs->create(nfr);
	fem.AddSurfaceLoad(pfs);
	
	// read the fluid flux data
	++tag;
	int nf[4], N;
	double s;
	for (int i=0; i<nfr; ++i)
	{
		FEFluidFlux::LOAD& fc = pfs->FluidFlux(i);
		FESurfaceElement& el = psurf->Element(i);
		
		sz = tag.AttributeValue("lc");
		fc.lc = atoi(sz) - 1;
		
		s  = atof(tag.AttributeValue("scale"));
		fc.s[0] = fc.s[1] = fc.s[2] = fc.s[3] = s;
		fc.s[4] = fc.s[5] = fc.s[6] = fc.s[7] = s;
		
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

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(pfs);
		pfs->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseBCSoluteFlux(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	
	const char* sz;
	bool blinear = false;
	sz = tag.AttributeValue("type", true);
	if (sz)
	{
		if (strcmp(sz, "linear") == 0) blinear = true;
		else if (strcmp(sz, "nonlinear") == 0) blinear = false;
		else throw XMLReader::InvalidAttributeValue(tag, "type", sz);
	}
	
	// count how many fluid flux cards there are
	int nfr = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(nfr);
	fem.GetMesh().AddSurface(psurf);
	
	// allocate fluid flux data
	FESoluteFlux* pfs = new FESoluteFlux(psurf, blinear);
	pfs->create(nfr);
	fem.AddSurfaceLoad(pfs);
	
	// read the fluid flux data
	++tag;
	int nf[4], N;
	double s;
	for (int i=0; i<nfr; ++i)
	{
		FESoluteFlux::LOAD& fc = pfs->SoluteFlux(i);
		FESurfaceElement& el = psurf->Element(i);
		
		sz = tag.AttributeValue("lc");
		fc.lc = atoi(sz) - 1;
		
		s  = atof(tag.AttributeValue("scale"));
		fc.s[0] = fc.s[1] = fc.s[2] = fc.s[3] = s;
		fc.s[4] = fc.s[5] = fc.s[6] = fc.s[7] = s;
		
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

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(pfs);
		pfs->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseBCHeatFlux(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// count how many heatflux cards there are
	int npr = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(npr);
	fem.GetMesh().AddSurface(psurf);

	// allocate flux data
	FEHeatFlux* ph = new FEHeatFlux(psurf);
	ph->create(npr);
	fem.AddSurfaceLoad(ph);

	const char* sz;

	// read the flux data
	++tag;
	int nf[4], N;
	double s;
	for (int i=0; i<npr; ++i)
	{
		FEHeatFlux::LOAD& pc = ph->HeatFlux(i);
		FESurfaceElement& el = psurf->Element(i);

		sz = tag.AttributeValue("lc");
		if (sz) pc.lc = atoi(sz) - 1;

		s  = atof(tag.AttributeValue("scale"));
		pc.s[0] = pc.s[1] = pc.s[2] = pc.s[3] = s;
		pc.s[4] = pc.s[5] = pc.s[6] = pc.s[7] = s;

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

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(ph);
		ph->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseBCConvectiveHeatFlux(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// count how many heatflux cards there are
	int npr = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(npr);
	fem.GetMesh().AddSurface(psurf);

	// allocate flux data
	FEConvectiveHeatFlux* ph = new FEConvectiveHeatFlux(psurf);
	ph->create(npr);
	fem.AddSurfaceLoad(ph);

	const char* sz;

	// read the flux data
	++tag;
	int nf[8], N;
	for (int i=0; i<npr; ++i)
	{
		FEConvectiveHeatFlux::LOAD& pc = ph->HeatFlux(i);
		FESurfaceElement& el = psurf->Element(i);

		sz = tag.AttributeValue("lc");
		if (sz) pc.lc = atoi(sz) - 1;

		double s  = atof(tag.AttributeValue("scale"));
		pc.s[0] = pc.s[1] = pc.s[2] = pc.s[3] = s;
		pc.s[4] = pc.s[5] = pc.s[6] = pc.s[7] = s;

		// heat transfer coefficient
		double hc = atof(tag.AttributeValue("hc"));
		pc.hc = hc;

		// read the element
		if      (tag == "quad4") el.SetType(FE_QUAD4G4);
		else if (tag == "tri3" ) el.SetType(m_pim->m_ntri3);
//		else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
//		else if (tag == "quad8") el.SetType(FE_QUAD8G9);
		else throw XMLReader::InvalidTag(tag);

		N = el.Nodes();
		tag.value(nf, N);
		for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;

		++tag;
	}

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(ph);
		ph->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseSpringSection(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// determine the spring type
	FEDiscreteMaterial* pm = 0;
	const char* szt = tag.AttributeValue("type", true);
	if (szt)
	{
		if (strcmp(szt, "linear") == 0) pm = new FELinearSpring;
		else if (strcmp(szt, "tension-only linear") == 0) pm = new FETensionOnlyLinearSpring;
		else if (strcmp(szt, "nonlinear") == 0) pm = new FENonLinearSpring;
	}
	else pm = new FELinearSpring;

	// create a new spring "domain"
	FEDiscreteSpringDomain* pd = new FEDiscreteSpringDomain(&mesh, pm);
	mesh.AddDomain(pd);

	pd->create(1);
	FEDiscreteElement& de = dynamic_cast<FEDiscreteElement&>(pd->ElementRef(0));
	de.SetType(FE_DISCRETE);
	de.m_nID = ++m_pim->m_maxid;
	
	// add a new material for each spring
	fem.AddMaterial(pm);
	pm->SetID(fem.Materials());
	de.SetMatID(fem.Materials()-1);

	int n[2];

	// read spring discrete elements
	++tag;
	do
	{
		if (tag == "node")
		{
			tag.value(n, 2);
			de.m_node[0] = n[0]-1;
			de.m_node[1] = n[1]-1;
		}
		else if (tag == "E") 
		{
			if (dynamic_cast<FELinearSpring*>(pm)) tag.value((dynamic_cast<FELinearSpring*>(pm))->m_E);
			else if (dynamic_cast<FETensionOnlyLinearSpring*>(pm)) tag.value((dynamic_cast<FETensionOnlyLinearSpring*>(pm))->m_E);
			else throw XMLReader::InvalidTag(tag);
		}
		else if (tag == "force")
		{
			if (dynamic_cast<FENonLinearSpring*>(pm))
			{
				FENonLinearSpring* ps = dynamic_cast<FENonLinearSpring*>(pm);
				tag.value(ps->m_F);
				const char* szl = tag.AttributeValue("lc");
				ps->m_nlc = atoi(szl) - 1;
			}
			else throw XMLReader::InvalidTag(tag);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	pd->InitMaterialPointData();
}

//-----------------------------------------------------------------------------
//! Parse the linear constraints section of the xml input file
//! This section is a subsection of the Boundary section

void FEBioBoundarySection::ParseConstraints(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// make sure there is a constraint defined
	if (tag.isleaf()) return;

	// read the master node
	FELinearConstraint LC;
	int node;
	tag.AttributeValue("node", node);
	LC.master.node = node-1;

	const char* szbc = tag.AttributeValue("bc");
	if      (strcmp(szbc, "x") == 0) LC.master.bc = 0;
	else if (strcmp(szbc, "y") == 0) LC.master.bc = 1;
	else if (strcmp(szbc, "z") == 0) LC.master.bc = 2;
	else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

	// we must deactive the master dof
	// so that it does not get assigned an equation
	fem.GetMesh().Node(node-1).m_ID[LC.master.bc] = -1;

	// read the slave nodes
	++tag;
	do
	{
		FELinearConstraint::SlaveDOF dof;
		if (tag == "node")
		{
			tag.value(dof.val);
			tag.AttributeValue("id", node);
			dof.node = node - 1;

			const char* szbc = tag.AttributeValue("bc");
			if      (strcmp(szbc, "x") == 0) dof.bc = 0;
			else if (strcmp(szbc, "y") == 0) dof.bc = 1;
			else if (strcmp(szbc, "z") == 0) dof.bc = 2;
			else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

			LC.slave.push_back(dof);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	// add the linear constraint to the system
	fem.m_LinC.push_back(LC);
}


//-----------------------------------------------------------------------------
//! Parses the contact section of the xml input file
//! The contact section is a subsection of the boundary section

void FEBioBoundarySection::ParseContactSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	// make sure that the version is 1.x
	int nversion = m_pim->Version();
	if (nversion >= 0x0200) throw XMLReader::InvalidTag(tag);

	const char* szt = tag.AttributeValue("type");

	if (strcmp(szt, "sliding_with_gaps") == 0)
	{
		// --- S L I D I N G   W I T H   G A P S ---

		FESlidingInterface* ps = new FESlidingInterface(&fem);
		fem.AddSurfacePairInteraction(ps);

		FEParameterList& pl = ps->GetParameterList();

		++tag;
		do
		{
			// read parameters
			if (m_pim->ReadParameter(tag, pl) == false)
			{
				if (tag == "surface")
				{
					const char* sztype = tag.AttributeValue("type");
					int ntype;
					if (strcmp(sztype, "master") == 0) ntype = 1;
					else if (strcmp(sztype, "slave") == 0) ntype = 2;

					FESlidingSurface& s = (ntype == 1? ps->m_ms : ps->m_ss);
					m.AddSurface(&s);

					int nfmt = 0;
					const char* szfmt = tag.AttributeValue("format", true);
					if (szfmt)
					{
						if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
						else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
					}

					// read the surface section
					ParseSurfaceSection(tag, s, nfmt, true);
				}
				else throw XMLReader::InvalidTag(tag);
			}

			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "facet-to-facet sliding") == 0)
	{
		// --- F A C E T   T O   F A C E T   S L I D I N G ---

		FEFacet2FacetSliding* ps = new FEFacet2FacetSliding(&fem);
		fem.AddSurfacePairInteraction(ps);

		FEParameterList& pl = ps->GetParameterList();

		++tag;
		do
		{
			// read parameters
			if (m_pim->ReadParameter(tag, pl) == false)
			{
				if (tag == "surface")
				{
					const char* sztype = tag.AttributeValue("type");
					int ntype;
					if (strcmp(sztype, "master") == 0) ntype = 1;
					else if (strcmp(sztype, "slave") == 0) ntype = 2;

					FEFacetSlidingSurface& s = (ntype == 1? ps->m_ms : ps->m_ss);
					m.AddSurface(&s);

					int nfmt = 0;
					const char* szfmt = tag.AttributeValue("format", true);
					if (szfmt)
					{
						if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
						else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
					}

					// read the surface section
					ParseSurfaceSection(tag, s, nfmt, false);
				}
				else throw XMLReader::InvalidTag(tag);
			}

			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "sliding-tension-compression") == 0)
	{
		// --- S L I D I N G   I N T E R F A C E   B W ---
		FESlidingInterfaceBW* ps = new FESlidingInterfaceBW(&fem);
		fem.AddSurfacePairInteraction(ps);
		
		FEParameterList& pl = ps->GetParameterList();
		
		++tag;
		do
		{
			// read parameters
			if (m_pim->ReadParameter(tag, pl) == false)
			{
				if (tag == "surface")
				{
					const char* sztype = tag.AttributeValue("type");
					int ntype;
					if (strcmp(sztype, "master") == 0) ntype = 1;
					else if (strcmp(sztype, "slave") == 0) ntype = 2;
					
					FESlidingSurfaceBW& s = (ntype == 1? ps->m_ms : ps->m_ss);
					m.AddSurface(&s);
					
					int nfmt = 0;
					const char* szfmt = tag.AttributeValue("format", true);
					if (szfmt)
					{
						if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
						else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
					}
					
					// read the surface section
					ParseSurfaceSection(tag, s, nfmt, false);
				}
				else throw XMLReader::InvalidTag(tag);
			}
			
			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "sliding2") == 0)
	{
		// --- S L I D I N G   I N T E R F A C E   2 ---
		FESlidingInterface2* ps = new FESlidingInterface2(&fem);
		fem.AddSurfacePairInteraction(ps);

		FEParameterList& pl = ps->GetParameterList();

		++tag;
		do
		{
			// read parameters
			if (m_pim->ReadParameter(tag, pl) == false)
			{
				if (tag == "surface")
				{
					const char* sztype = tag.AttributeValue("type");
					int ntype;
					if (strcmp(sztype, "master") == 0) ntype = 1;
					else if (strcmp(sztype, "slave") == 0) ntype = 2;

					FESlidingSurface2& s = (ntype == 1? ps->m_ms : ps->m_ss);
					m.AddSurface(&s);

					int nfmt = 0;
					const char* szfmt = tag.AttributeValue("format", true);
					if (szfmt)
					{
						if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
						else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
					}

					// read the surface section
					ParseSurfaceSection(tag, s, nfmt, false);
				}
				else throw XMLReader::InvalidTag(tag);
			}

			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "sliding3") == 0)
	{
		// --- S L I D I N G   I N T E R F A C E   3 ---
		FESlidingInterface3* ps = new FESlidingInterface3(&fem);
		fem.AddSurfacePairInteraction(ps);

		FEParameterList& pl = ps->GetParameterList();
		
		++tag;
		do
		{
			// read parameters
			if (m_pim->ReadParameter(tag, pl) == false)
			{
				if (tag == "surface")
				{
					const char* sztype = tag.AttributeValue("type");
					int ntype;
					if (strcmp(sztype, "master") == 0) ntype = 1;
					else if (strcmp(sztype, "slave") == 0) ntype = 2;
					
					FESlidingSurface3& s = (ntype == 1? ps->m_ms : ps->m_ss);
					m.AddSurface(&s);
					
					int nfmt = 0;
					const char* szfmt = tag.AttributeValue("format", true);
					if (szfmt)
					{
						if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
						else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
					}
					
					// read the surface section
					ParseSurfaceSection(tag, s, nfmt, false);
				}
				else throw XMLReader::InvalidTag(tag);
			}
			
			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "tied") == 0)
	{
		// --- T I E D   C O N T A C T  ---

		FETiedInterface* ps = new FETiedInterface(&fem);
		fem.AddSurfacePairInteraction(ps);

		// add this contact interface to the current step
		if (m_pim->m_nsteps > 0)
		{
			GetStep()->AddSurfacePairInteraction(ps);
			ps->Deactivate();
		}

		FEParameterList& pl = ps->GetParameterList();

		++tag;
		do
		{
			if (m_pim->ReadParameter(tag, pl) == false)
			{
				if (tag == "surface")
				{
					const char* sztype = tag.AttributeValue("type");
					int ntype;
					if (strcmp(sztype, "master") == 0) ntype = 1;
					else if (strcmp(sztype, "slave") == 0) ntype = 2;

					FETiedContactSurface& s = (ntype == 1? ps->ms : ps->ss);
					m.AddSurface(&s);

					int nfmt = 0;
					const char* szfmt = tag.AttributeValue("format", true);
					if (szfmt)
					{
						if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
						else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
					}

					// read the surface section
					ParseSurfaceSection(tag, s, nfmt, true);
				}
				else throw XMLReader::InvalidTag(tag);
			}

			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "facet-to-facet tied") == 0)
	{
		// --- F2F T I E D   C O N T A C T  ---

		FEFacet2FacetTied* pt = new FEFacet2FacetTied(&fem);
		fem.AddSurfacePairInteraction(pt);

		// add this contact interface to the current step
		if (m_pim->m_nsteps > 0)
		{
			GetStep()->AddSurfacePairInteraction(pt);
			pt->Deactivate();
		}

		FEParameterList& pl = pt->GetParameterList();

		++tag;
		do
		{
			if (m_pim->ReadParameter(tag, pl) == false)
			{
				if (tag == "surface")
				{
					const char* sztype = tag.AttributeValue("type");
					int ntype;
					if (strcmp(sztype, "master") == 0) ntype = 1;
					else if (strcmp(sztype, "slave") == 0) ntype = 2;

					FEFacetTiedSurface& s = dynamic_cast<FEFacetTiedSurface&>((ntype == 1? *pt->GetMasterSurface(): *pt->GetSlaveSurface()));
					m.AddSurface(&s);

					int nfmt = 0;
					const char* szfmt = tag.AttributeValue("format", true);
					if (szfmt)
					{
						if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
						else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
					}

					// read the surface section
					ParseSurfaceSection(tag, s, nfmt, true);
				}
				else throw XMLReader::InvalidTag(tag);
			}

			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "tied-biphasic") == 0)
	{
		// --- T I E D - B I P H A S I C ---
		FETiedBiphasicInterface* ps = new FETiedBiphasicInterface(&fem);
		fem.AddSurfacePairInteraction(ps);
		
		FEParameterList& pl = ps->GetParameterList();
		
		++tag;
		do
		{
			// read parameters
			if (m_pim->ReadParameter(tag, pl) == false)
			{
				if (tag == "surface")
				{
					const char* sztype = tag.AttributeValue("type");
					int ntype;
					if (strcmp(sztype, "master") == 0) ntype = 1;
					else if (strcmp(sztype, "slave") == 0) ntype = 2;
					
					FETiedBiphasicSurface& s = (ntype == 1? ps->m_ms : ps->m_ss);
					m.AddSurface(&s);
					
					int nfmt = 0;
					const char* szfmt = tag.AttributeValue("format", true);
					if (szfmt)
					{
						if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
						else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
					}
					
					// read the surface section
					ParseSurfaceSection(tag, s, nfmt, false);
				}
				else throw XMLReader::InvalidTag(tag);
			}
			
			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "periodic boundary") == 0)
	{
		// --- P E R I O D I C   B O U N D A R Y  ---

		FEPeriodicBoundary* ps = new FEPeriodicBoundary(&fem);
		fem.AddSurfacePairInteraction(ps);

		FEParameterList& pl = ps->GetParameterList();

		++tag;
		do
		{
			if (m_pim->ReadParameter(tag, pl) == false)
			{
				if (tag == "surface")
				{
					const char* sztype = tag.AttributeValue("type");
					int ntype;
					if (strcmp(sztype, "master") == 0) ntype = 1;
					else if (strcmp(sztype, "slave") == 0) ntype = 2;

					FEPeriodicSurface& s = (ntype == 1? ps->m_ms : ps->m_ss);
					m.AddSurface(&s);

					int nfmt = 0;
					const char* szfmt = tag.AttributeValue("format", true);
					if (szfmt)
					{
						if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
						else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
					}

					// read the surface section
					ParseSurfaceSection(tag, s, nfmt, true);
				}
				else throw XMLReader::InvalidTag(tag);
			}

			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "surface constraint") == 0)
	{
		// --- S U R F A C E   C O N S T R A I N T ---

		FESurfaceConstraint* ps = new FESurfaceConstraint(&fem);
		fem.AddSurfacePairInteraction(ps);

		FEParameterList& pl = ps->GetParameterList();

		++tag;
		do
		{
			if (m_pim->ReadParameter(tag, pl) == false)
			{
				if (tag == "surface")
				{
					const char* sztype = tag.AttributeValue("type");
					int ntype;
					if (strcmp(sztype, "master") == 0) ntype = 1;
					else if (strcmp(sztype, "slave") == 0) ntype = 2;

					FESurfaceConstraintSurface& s = (ntype == 1? ps->m_ms : ps->m_ss);
					m.AddSurface(&s);

					int nfmt = 0;
					const char* szfmt = tag.AttributeValue("format", true);
					if (szfmt)
					{
						if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
						else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
					}

					// read the surface section
					ParseSurfaceSection(tag, s, nfmt, true);
				}
				else throw XMLReader::InvalidTag(tag);
			}

			++tag;
		}
		while (!tag.isend());
	}	
	else if (strcmp(szt, "rigid_wall") == 0)
	{
		// --- R I G I D   W A L L   I N T E R F A C E ---

		FERigidWallInterface* ps = new FERigidWallInterface(&fem);
		fem.AddSurfacePairInteraction(ps);

		FEParameterList& pl = ps->GetParameterList();

		++tag;
		do
		{
			if (m_pim->ReadParameter(tag, pl) == false)
			{
				if (tag == "plane")
				{
					ps->SetMasterSurface(new FEPlane(&fem));
					FEPlane& pl = dynamic_cast<FEPlane&>(*ps->m_mp);
					const char* sz = tag.AttributeValue("lc", true);
					if (sz)	pl.m_nplc = atoi(sz) - 1;

					double* a = pl.GetEquation();
					tag.value(a, 4);
				}
				else if (tag == "sphere")
				{
					ps->SetMasterSurface(new FERigidSphere(&fem));
					FERigidSphere& s = dynamic_cast<FERigidSphere&>(*ps->m_mp);
					++tag;
					do
					{
						if      (tag == "center") tag.value(s.m_rc);
						else if (tag == "radius") tag.value(s.m_R);
						else if (tag == "xtrans")
						{
							const char* szlc = tag.AttributeValue("lc");
							s.m_nplc[0] = atoi(szlc) - 1;
						}
						else if (tag == "ytrans")
						{
							const char* szlc = tag.AttributeValue("lc");
							s.m_nplc[1] = atoi(szlc) - 1;
						}
						else if (tag == "ztrans")
						{
							const char* szlc = tag.AttributeValue("lc");
							s.m_nplc[2] = atoi(szlc) - 1;
						}
						else throw XMLReader::InvalidTag(tag);
						++tag;
					}
					while (!tag.isend());
				}
				else if (tag == "surface")
				{
					FERigidWallSurface& s = ps->m_ss;

					int nfmt = 0;
					const char* szfmt = tag.AttributeValue("format", true);
					if (szfmt)
					{
						if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
						else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
					}

					// read the surface section
					ParseSurfaceSection(tag, s, nfmt, true);
				}
				else throw XMLReader::InvalidTag(tag);
			}
			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "rigid") == 0)
	{
		// --- R I G I D   B O D Y   I N T E R F A C E ---

		// count how many rigid nodes there are
		int nrn= 0;
		XMLTag t(tag); ++t;
		while (!t.isend()) { nrn++; ++t; }

		++tag;
		int id, rb;
		for (int i=0; i<nrn; ++i)
		{
			id = atoi(tag.AttributeValue("id"))-1;
			rb = atoi(tag.AttributeValue("rb"))-1;

			FERigidNode* prn = new FERigidNode;

			prn->nid = id;
			prn->rid = rb;
			fem.AddRigidNode(prn);

			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddBoundaryCondition(prn);
				prn->Deactivate();
			}

			++tag;
		}
	}
	else if (strcmp(szt, "rigid joint") == 0)
	{
		// --- R I G I D   J O I N T   I N T E R F A C E ---

		FERigidJoint* prj = new FERigidJoint(&fem);
		FEParameterList& pl = prj->GetParameterList();
		++tag;
		do
		{
			if (m_pim->ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
		prj->m_nRBa--;
		prj->m_nRBb--;
		fem.AddNonlinearConstraint(prj);
	}
	else if (strcmp(szt, "linear constraint") == 0)
	{
		FEModel& fem = *GetFEModel();

		// make sure there is a constraint defined
		if (tag.isleaf()) return;

		// create a new linear constraint manager
		FELinearConstraintSet* pLCS = new FELinearConstraintSet(&fem);
		fem.AddNonlinearConstraint(pLCS);

		// read the linear constraints
		++tag;
		do
		{
			if (tag == "linear_constraint")
			{
				FEAugLagLinearConstraint* pLC = new FEAugLagLinearConstraint;

				FEAugLagLinearConstraint::DOF dof;
				++tag;
				do
				{
					if (tag == "node")
					{
						tag.value(dof.val);
						int node;
						tag.AttributeValue("id", node);
						dof.node = node - 1;

						const char* szbc = tag.AttributeValue("bc");
						if      (strcmp(szbc, "x") == 0) dof.bc = 0;
						else if (strcmp(szbc, "y") == 0) dof.bc = 1;
						else if (strcmp(szbc, "z") == 0) dof.bc = 2;
						else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

						pLC->m_dof.push_back(dof);
					}
					else throw XMLReader::InvalidTag(tag);
					++tag;
				}
				while (!tag.isend());

				// add the linear constraint to the system
				pLCS->add(pLC);
			}
			else if (tag == "tol"    ) tag.value(pLCS->m_tol);
			else if (tag == "penalty") tag.value(pLCS->m_eps);
			else if (tag == "maxaug") tag.value(pLCS->m_naugmax);
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}
	else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
}
