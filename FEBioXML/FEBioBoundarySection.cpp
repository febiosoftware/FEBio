#include "stdafx.h"
#include "FEBioBoundarySection.h"
#include "FECore/FEDiscreteMaterial.h"
#include "FECore/FEDiscreteDomain.h"
#include "FEBioMech/FERigidWallInterface.h"
#include "FEBioMech/FEAugLagLinearConstraint.h"
#include "FEBioMech/FERigidJoint.h"
#include "FEBioMech/FERigidSphericalJoint.h"
#include "FEBioMech/FERigidRevoluteJoint.h"
#include "FEBioMech/FERigidPrismaticJoint.h"
#include "FEBioMech/FERigidCylindricalJoint.h"
#include "FEBioMech/FERigidPlanarJoint.h"
#include "FEBioMech/FERigidSpring.h"
#include "FEBioMech/FERigidDamper.h"
#include "FEBioMech/FERigidAngularDamper.h"
#include "FEBioMech/FERigidContractileForce.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
//!  Parses the boundary section from the xml file
//!
void FEBioBoundarySection::Parse(XMLTag& tag)
{
	// make sure this tag has children (for older versions)
	if (m_pim->Version() < 0x0200) 
		if (tag.isleaf()) return;

	++tag;
	do
	{
		if (m_pim->Version() < 0x0200)
		{
			if      (tag == "fix"              ) ParseBCFix         (tag);
			else if (tag == "prescribe"        ) ParseBCPrescribe   (tag);
			else if (tag == "contact"          ) ParseContactSection(tag);
			else if (tag == "linear_constraint") ParseConstraints   (tag);
			else if (tag == "spring"           ) ParseSpringSection (tag);
			else throw XMLReader::InvalidTag(tag);
		}
		else
		{
			if      (tag == "fix"              ) ParseBCFix20      (tag);
			else if (tag == "prescribe"        ) ParseBCPrescribe20(tag);
			else if (tag == "linear_constraint") ParseConstraints  (tag);
			else throw XMLReader::InvalidTag(tag);
		}
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

	int N, nf[9];

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
			else if (tag == "tri6" ) el.SetType(FE_TRI6NI );
			else throw XMLReader::InvalidTag(tag);
		}
		else
		{
			if      (tag == "quad4") el.SetType(FE_QUAD4G4);
			else if (tag == "tri3" ) el.SetType(m_pim->m_ntri3);
			else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
			else if (tag == "tri7" ) el.SetType(m_pim->m_ntri7);
			else if (tag == "quad8") el.SetType(FE_QUAD8G9);
			else if (tag == "quad9") el.SetType(FE_QUAD9G9);
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
				int ne[9];
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
			if      (strcmp(sz, "x"  ) == 0) { node.m_BC[DOF_X] = -1; }
			else if (strcmp(sz, "y"  ) == 0) { node.m_BC[DOF_Y] = -1; }
			else if (strcmp(sz, "z"  ) == 0) { node.m_BC[DOF_Z] = -1; }
			else if (strcmp(sz, "xy" ) == 0) { node.m_BC[DOF_X] = node.m_BC[DOF_Y] = -1; }
			else if (strcmp(sz, "yz" ) == 0) { node.m_BC[DOF_Y] = node.m_BC[DOF_Z] = -1; }
			else if (strcmp(sz, "xz" ) == 0) { node.m_BC[DOF_X] = node.m_BC[DOF_Z] = -1; }
			else if (strcmp(sz, "xyz") == 0) { node.m_BC[DOF_X] = node.m_BC[DOF_Y] = node.m_BC[DOF_Z] = -1; }
			else if (strcmp(sz, "p"  ) == 0) { node.m_BC[DOF_P] = -1; }
			else if (strcmp(sz, "u"  ) == 0) { node.m_BC[DOF_U] = -1; }
			else if (strcmp(sz, "v"  ) == 0) { node.m_BC[DOF_V] = -1; }
			else if (strcmp(sz, "w"  ) == 0) { node.m_BC[DOF_W] = -1; }
			else if (strcmp(sz, "uv" ) == 0) { node.m_BC[DOF_U] = node.m_BC[DOF_V] = -1; }
			else if (strcmp(sz, "vw" ) == 0) { node.m_BC[DOF_V] = node.m_BC[DOF_W] = -1; }
			else if (strcmp(sz, "uw" ) == 0) { node.m_BC[DOF_U] = node.m_BC[DOF_W] = -1; }
			else if (strcmp(sz, "uvw") == 0) { node.m_BC[DOF_U] = node.m_BC[DOF_V] = node.m_BC[DOF_W] = -1; }
			else if (strcmp(sz, "t"  ) == 0) { node.m_BC[DOF_T] = -1; }
			else if (strcmp(sz, "c"  ) == 0) { node.m_BC[DOF_C] = -1; }
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
			if      (strcmp(sz, "x"  ) == 0) { node.m_BC[DOF_X] = -1; }
			else if (strcmp(sz, "y"  ) == 0) { node.m_BC[DOF_Y] = -1; }
			else if (strcmp(sz, "z"  ) == 0) { node.m_BC[DOF_Z] = -1; }
			else if (strcmp(sz, "xy" ) == 0) { node.m_BC[DOF_X] = node.m_BC[DOF_Y] = -1; }
			else if (strcmp(sz, "yz" ) == 0) { node.m_BC[DOF_Y] = node.m_BC[DOF_Z] = -1; }
			else if (strcmp(sz, "xz" ) == 0) { node.m_BC[DOF_X] = node.m_BC[DOF_Z] = -1; }
			else if (strcmp(sz, "xyz") == 0) { node.m_BC[DOF_X] = node.m_BC[DOF_Y] = node.m_BC[DOF_Z] = -1; }
			else if (strcmp(sz, "p"  ) == 0) { node.m_BC[DOF_P] = -1; }
			else if (strcmp(sz, "u"  ) == 0) { node.m_BC[DOF_U] = -1; }
			else if (strcmp(sz, "v"  ) == 0) { node.m_BC[DOF_V] = -1; }
			else if (strcmp(sz, "w"  ) == 0) { node.m_BC[DOF_W] = -1; }
			else if (strcmp(sz, "uv" ) == 0) { node.m_BC[DOF_U] = node.m_BC[DOF_V] = -1; }
			else if (strcmp(sz, "vw" ) == 0) { node.m_BC[DOF_V] = node.m_BC[DOF_W] = -1; }
			else if (strcmp(sz, "uw" ) == 0) { node.m_BC[DOF_U] = node.m_BC[DOF_W] = -1; }
			else if (strcmp(sz, "uvw") == 0) { node.m_BC[DOF_U] = node.m_BC[DOF_V] = node.m_BC[DOF_W] = -1; }
			else if (strcmp(sz, "t"  ) == 0) { node.m_BC[DOF_T] = -1; }
			else if (strcmp(sz, "c"  ) == 0) { node.m_BC[DOF_C] = -1; }
			else if (strncmp(sz, "c", 1) == 0) node.m_BC[DOF_C + atoi(&sz[1]) - 1] = -1;
			else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);
			++tag;
		}
		while (!tag.isend());
	}
}


//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseBCFix20(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NN = mesh.Nodes();

	// make sure this section does not appear in a step section
	if (m_pim->m_nsteps != 0) throw XMLReader::InvalidTag(tag);

	// get the required bc attribute
	char szbc[8];
	strcpy(szbc, tag.AttributeValue("bc"));

	// process the bc string
	vector<int> bc;
	char* ch = szbc;
	while (*ch)
	{
		if      (*ch == 'x') bc.push_back(DOF_X);
		else if (*ch == 'y') bc.push_back(DOF_Y);
		else if (*ch == 'z') bc.push_back(DOF_Z);
		else if (*ch == 'u') bc.push_back(DOF_U);
		else if (*ch == 'v') bc.push_back(DOF_V);
		else if (*ch == 'w') bc.push_back(DOF_W);
		else if (*ch == 'p') bc.push_back(DOF_P);
		else if (*ch == 't') bc.push_back(DOF_T);
		else if (*ch == 'c')
		{
			char ci = ch[1];
			if (isdigit(ci))
			{
				int i = (ci - '0');
				bc.push_back(DOF_C + i - 1);
				++ch;
			}
			else bc.push_back(DOF_C);
		}
		else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);
		if (*ch) ++ch;
	}
	if (bc.empty()) throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);
	int nbc = bc.size();

	// see if the set attribute is defined
	const char* szset = tag.AttributeValue("set", true);
	if (szset)
	{
		// make sure the tag is a leaf
		if (tag.isleaf() == false) throw XMLReader::InvalidValue(tag);

		// process the node set
		FENodeSet* pns = mesh.FindNodeSet(szset);
		if (pns == 0) throw XMLReader::InvalidAttributeValue(tag, "set", szset);

		FENodeSet& ns = *pns;
		int N = pns->size();
		for (int i=0; i<N; ++i)
		{
			int n = ns[i];
			if ((n < 0)||(n >= NN)) throw XMLReader::InvalidTag(tag);
			FENode& node = mesh.Node(n);
			for (int j=0; j<nbc; ++j) node.m_BC[bc[j]] = -1;
		}
	}
	else
	{
		// Read the fixed nodes
		++tag;
		do
		{
			int n = atoi(tag.AttributeValue("id"))-1;
			if ((n<0) || (n >= NN)) throw XMLReader::InvalidAttributeValue(tag, "id");
			FENode& node = fem.GetMesh().Node(n);
			for (int j=0; j<nbc; ++j) node.m_BC[bc[j]] = -1;
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

	assert(nversion < 0x0200);

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
		m_pim->value(tag, s);

		// loop over all nodes in the nodeset
		FENodeSet& ns = *ps;
		int N = ns.size();
		for (int i=0; i<N; ++i)
		{
			FEPrescribedBC* pdc = new FEPrescribedBC(&fem);
			pdc->node = ns[i];
			pdc->bc = bc;
			pdc->lc = lc;
			pdc->s = s;
			fem.AddPrescribedBC(pdc);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddModelComponent(pdc);
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

			FEPrescribedBC* pdc = new FEPrescribedBC(&fem);
			pdc->node = n;
			pdc->bc = bc;
			pdc->lc = lc;
			tag.value(pdc->s);
			pdc->br = br;
			fem.AddPrescribedBC(pdc);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddModelComponent(pdc);
				pdc->Deactivate();
			}
			++tag;
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseBCPrescribe20(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NN = mesh.Nodes();

	int nversion = m_pim->Version();

	assert(nversion >= 0x0200);

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

	// get the load curve number
	sz = tag.AttributeValue("lc");
	int lc = atoi(sz) - 1;

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
		// TODO: the scale attribute is obsolete. Use the tag's value instead
		// TODO: instead of defining everything in attributes, perhaps I should
		//       just make child elements.
		double scale = 1.0;
		tag.AttributeValue("scale", scale, true);
		if (tag.isempty() == false) m_pim->value(tag, scale);

		FENodeSet& ns = *pns;
		int N = ns.size();
		for (int i=0; i<N; ++i)
		{
			int n = ns[i];
			FEPrescribedBC* pbc = new FEPrescribedBC(&fem);
			pbc->node = n;
			pbc->bc = bc;
			pbc->lc = lc;
			pbc->s = scale;
			pbc->br = br;
			fem.AddPrescribedBC(pbc);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddModelComponent(pbc);
				pbc->Deactivate();
			}
		}
	}
	else
	{
		// read the prescribed data
		++tag;
		for (int i=0; i<ndis; ++i)
		{
			// get the node ID
			int n = atoi(tag.AttributeValue("id"))-1;

			// create a new BC
			FEPrescribedBC* pdc = new FEPrescribedBC(&fem);
			pdc->node = n;
			pdc->bc = bc;
			pdc->lc = lc;
			pdc->br = br;
			m_pim->value(tag, pdc->s);
			fem.AddPrescribedBC(pdc);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddModelComponent(pdc);
				pdc->Deactivate();
			}
			++tag;
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseSpringSection(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// determine the spring type
	const char* szt = tag.AttributeValue("type", true);
	if (szt == 0) szt = "linear";
	FEDiscreteMaterial* pm = dynamic_cast<FEDiscreteMaterial*>(fecore_new<FEMaterial>(FEMATERIAL_ID, szt, &fem));
	if (pm == 0) throw XMLReader::InvalidAttributeValue(tag, "type", szt);

	// create a new spring "domain"
	FECoreKernel& febio = FECoreKernel::GetInstance();
	FE_Element_Spec spec;
	spec.eshape = ET_TRUSS2;
	spec.etype  = FE_DISCRETE;
	int ndomtype = febio.GetDomainType(spec, pm);
	FEDiscreteDomain* pd = dynamic_cast<FEDiscreteDomain*>(febio.CreateDomain(ndomtype, &mesh, pm));
	mesh.AddDomain(pd);

	pd->create(1);
	FEDiscreteElement& de = pd->Element(0);
	de.SetType(FE_DISCRETE);
	de.m_nID = ++m_pim->m_maxid;
	
	// add a new material for each spring
	fem.AddMaterial(pm);
	pm->SetID(fem.Materials());
	de.SetMatID(fem.Materials()-1);

	// read spring discrete elements
	++tag;
	do
	{
		// read the required node tag
		if (tag == "node")
		{
			int n[2];
			tag.value(n, 2);
			de.m_node[0] = n[0]-1;
			de.m_node[1] = n[1]-1;
		}
		else
		{
			// read the actual spring material parameters
			FEParameterList& pl = pm->GetParameterList();
			if (m_pim->ReadParameter(tag, pl) == 0)
			{
				throw XMLReader::InvalidTag(tag);
			}
		}
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
	fem.GetMesh().Node(node-1).m_BC[LC.master.bc] = -1;

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
void FEBioBoundarySection::ParseContactInterface(XMLTag& tag, FESurfacePairInteraction* pci)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	// get the parameter list
	FEParameterList& pl = pci->GetParameterList();

	// read the parameters
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

				FESurface& s = *(ntype == 1? pci->GetMasterSurface() : pci->GetSlaveSurface());
				m.AddSurface(&s);

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}

				// read the surface section
				ParseSurfaceSection(tag, s, nfmt, pci->UseNodalIntegration());
			}
			else throw XMLReader::InvalidTag(tag);
		}

		++tag;
	}
	while (!tag.isend());
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
/*
	// TODO: This was defined in the tied interface parse
	//		 Not sure what to do with this.

		// add this contact interface to the current step
		if (m_pim->m_nsteps > 0)
		{
			GetStep()->AddSurfacePairInteraction(ps);
			ps->Deactivate();
		}
*/

	// get the type attribute
	const char* szt = tag.AttributeValue("type");

	// Not all contact interfaces can be parsed automatically.
	// First, check all these special cases.
	if (strcmp(szt, "rigid_wall") == 0)
	{
		// --- R I G I D   W A L L   I N T E R F A C E ---

		FERigidWallInterface* ps = dynamic_cast<FERigidWallInterface*>(fecore_new<FESurfacePairInteraction>(FESURFACEPAIRINTERACTION_ID, szt, GetFEModel()));
		if (ps)
		{
			fem.AddSurfacePairInteraction(ps);

			++tag;
			do
			{
				if (m_pim->ReadParameter(tag, ps) == false)
				{
					if (tag == "surface")
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
		else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
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

			FERigidNode* prn = new FERigidNode(&fem);

			prn->nid = id;
			prn->rid = rb;
			fem.AddRigidNode(prn);

			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddModelComponent(prn);
				prn->Deactivate();
			}

			++tag;
		}
	}
	else if (strcmp(szt, "rigid joint") == 0)
	{
		// --- R I G I D   J O I N T   I N T E R F A C E ---

		FERigidJoint* prj = dynamic_cast<FERigidJoint*>(fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, szt, GetFEModel()));
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
	else if (strcmp(szt, "rigid spherical joint") == 0)
	{
		// --- R I G I D   S P H E R I C A L   J O I N T   I N T E R F A C E ---
        
		FERigidSphericalJoint* prj = dynamic_cast<FERigidSphericalJoint*>(fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, szt, GetFEModel()));
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
    else if (strcmp(szt, "rigid revolute joint") == 0)
    {
        // --- R I G I D   R E V O L U T E  J O I N T   I N T E R F A C E ---
        
        FERigidRevoluteJoint* prj = dynamic_cast<FERigidRevoluteJoint*>(fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, szt, GetFEModel()));
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
    else if (strcmp(szt, "rigid prismatic joint") == 0)
    {
        // --- R I G I D   P R I S M A T I C  J O I N T   I N T E R F A C E ---
        
        FERigidPrismaticJoint* prj = dynamic_cast<FERigidPrismaticJoint*>(fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, szt, GetFEModel()));
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
    else if (strcmp(szt, "rigid cylindrical joint") == 0)
    {
        // --- R I G I D   C Y L I N D R I C A L  J O I N T   I N T E R F A C E ---
        
        FERigidCylindricalJoint* prj = dynamic_cast<FERigidCylindricalJoint*>(fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, szt, GetFEModel()));
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
    else if (strcmp(szt, "rigid planar joint") == 0)
    {
        // --- R I G I D   P L A N A R  J O I N T   I N T E R F A C E ---
        
        FERigidPlanarJoint* prj = dynamic_cast<FERigidPlanarJoint*>(fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, szt, GetFEModel()));
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
    else if (strcmp(szt, "rigid spring") == 0)
    {
        // --- R I G I D   S P R I N G   I N T E R F A C E ---
        
        FERigidSpring* prj = dynamic_cast<FERigidSpring*>(fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, szt, GetFEModel()));
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
    else if (strcmp(szt, "rigid damper") == 0)
    {
        // --- R I G I D   D A M P E R   I N T E R F A C E ---
        
        FERigidDamper* prj = dynamic_cast<FERigidDamper*>(fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, szt, GetFEModel()));
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
    else if (strcmp(szt, "rigid angular damper") == 0)
    {
        // --- R I G I D   A N G U L A R D A M P E R   I N T E R F A C E ---
        
        FERigidAngularDamper* prj = dynamic_cast<FERigidAngularDamper*>(fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, szt, GetFEModel()));
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
    else if (strcmp(szt, "rigid contractile force") == 0)
    {
        // --- R I G I D   C O N T R A C T I L E F O R C E   I N T E R F A C E ---
        
        FERigidContractileForce* prj = dynamic_cast<FERigidContractileForce*>(fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, szt, GetFEModel()));
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
		FELinearConstraintSet* pLCS = dynamic_cast<FELinearConstraintSet*>(fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, szt, GetFEModel()));
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
			else if (tag == "tol"    ) m_pim->value(tag, pLCS->m_tol);
			else if (tag == "penalty") m_pim->value(tag, pLCS->m_eps);
			else if (tag == "maxaug" ) m_pim->value(tag, pLCS->m_naugmax);
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}
	else
	{
		// If we get here, we try to create a contact interface
		// using the FEBio kernel. 
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fecore_new<FESurfacePairInteraction>(FESURFACEPAIRINTERACTION_ID, szt, GetFEModel()));
		if (pci)
		{
			fem.AddSurfacePairInteraction(pci);
			ParseContactInterface(tag, pci);
			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddModelComponent(pci);
				pci->Deactivate();
			}
		}
		else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
	}
}
