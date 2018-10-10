#include "stdafx.h"
#include "FEBioBoundarySection.h"
#include <FECore/FEModel.h>
#include <FECore/FEDiscreteMaterial.h>
#include <FECore/FEDiscreteDomain.h>
#include <FEBioMech/FERigidWallInterface.h>
#include <FEBioMech/FEAugLagLinearConstraint.h>
#include <FEBioMech/FERigidForce.h>
#include <FECore/FEPrescribedDOF.h>
#include <FECore/FEFixedBC.h>
#include <FEBioMech/RigidBC.h>
#include <FEBioMech/FERigidSystem.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FELinearConstraintManager.h>
#include <FEBioMech/FEPeriodicLinearConstraint.h>
#include <FEBioMech/FEPeriodicLinearConstraint2O.h>
#include <FEBioMech/FEMergedConstraint.h>
#include <FEBioMech/FEMechModel.h>
#include <FEBioMech/FERigidMaterial.h>

//-----------------------------------------------------------------------------
void FEBioBoundarySection1x::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "fix"              ) ParseBCFix         (tag);
		else if (tag == "prescribe"        ) ParseBCPrescribe   (tag);
		else if (tag == "contact"          ) ParseContactSection(tag);
		else if (tag == "linear_constraint") ParseConstraints   (tag);
		else if (tag == "spring"           ) ParseSpringSection (tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection2::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "fix"              ) ParseBCFix      (tag);
		else if (tag == "prescribe"        ) ParseBCPrescribe(tag);
		else if (tag == "linear_constraint") ParseConstraints(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection25::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	// build the node set map for faster lookup
	BuildNodeSetMap();

	++tag;
	do
	{
		if      (tag == "fix"                          ) ParseBCFix                     (tag);
		else if (tag == "prescribe"                    ) ParseBCPrescribe               (tag);
		else if (tag == "bc"                           ) ParseBC                        (tag);
		else if (tag == "rigid"                        ) ParseBCRigid                   (tag);
		else if (tag == "rigid_body"                   ) ParseRigidBody                 (tag);
		else if (tag == "linear_constraint"            ) ParseConstraints               (tag);
		else if (tag == "periodic_linear_constraint"   ) ParsePeriodicLinearConstraint  (tag);
		else if (tag == "periodic_linear_constraint_2O") ParsePeriodicLinearConstraint2O(tag);
		else if (tag == "merge"                        ) ParseMergeConstraint           (tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//---------------------------------------------------------------------------------
void FEBioBoundarySection25::BuildNodeSetMap()
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	m_NodeSet.clear();
	for (int i=0; i<m.NodeSets(); ++i)
	{
		FENodeSet* nsi = m.NodeSet(i);
		m_NodeSet[nsi->GetName()] = nsi;
	}
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
	s.Create(faces);

	FEModelBuilder* feb = GetBuilder();

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
            else if (tag == "quad8" ) el.SetType(FE_QUAD8NI);
            else if (tag == "quad9" ) el.SetType(FE_QUAD9NI);
			else throw XMLReader::InvalidTag(tag);
		}
		else
		{
			if      (tag == "quad4") el.SetType(FE_QUAD4G4);
			else if (tag == "tri3" ) el.SetType(feb->m_ntri3);
			else if (tag == "tri6" ) el.SetType(feb->m_ntri6);
			else if (tag == "tri7" ) el.SetType(feb->m_ntri7);
			else if (tag == "tri10") el.SetType(feb->m_ntri10);
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
				el.m_elem[0] = pe;
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
	DOFS& dofs = fem.GetDOFS();

	// get the DOF indices
	const int dof_X = fem.GetDOFIndex("x");
	const int dof_Y = fem.GetDOFIndex("y");
	const int dof_Z = fem.GetDOFIndex("z");
	const int dof_U = fem.GetDOFIndex("u");
	const int dof_V = fem.GetDOFIndex("v");
	const int dof_W = fem.GetDOFIndex("w");
    const int dof_SX = fem.GetDOFIndex("sx");
    const int dof_SY = fem.GetDOFIndex("sy");
    const int dof_SZ = fem.GetDOFIndex("sz");
    const int dof_WX = fem.GetDOFIndex("wx");
    const int dof_WY = fem.GetDOFIndex("wy");
    const int dof_WZ = fem.GetDOFIndex("wz");

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
			int n = s[i];

			int ndof = dofs.GetDOF(sz);
			if (ndof >= 0) fem.AddFixedBC(n, ndof);
			else
			{
				// The supported fixed BC strings don't quite follow the dof naming convention.
				// For now, we'll check these BC explicitly, but I want to get rid of this in the future.
				if      (strcmp(sz, "xy" ) == 0) { fem.AddFixedBC(n, dof_X); fem.AddFixedBC(n, dof_Y); }
				else if (strcmp(sz, "yz" ) == 0) { fem.AddFixedBC(n, dof_Y); fem.AddFixedBC(n, dof_Z); }
				else if (strcmp(sz, "xz" ) == 0) { fem.AddFixedBC(n, dof_X); fem.AddFixedBC(n, dof_Z); }
				else if (strcmp(sz, "xyz") == 0) { fem.AddFixedBC(n, dof_X); fem.AddFixedBC(n, dof_Y); fem.AddFixedBC(n, dof_Z); }
				else if (strcmp(sz, "uv" ) == 0) { fem.AddFixedBC(n, dof_U); fem.AddFixedBC(n, dof_V); }
				else if (strcmp(sz, "vw" ) == 0) { fem.AddFixedBC(n, dof_V); fem.AddFixedBC(n, dof_W); }
				else if (strcmp(sz, "uw" ) == 0) { fem.AddFixedBC(n, dof_U); fem.AddFixedBC(n, dof_W); }
				else if (strcmp(sz, "uvw") == 0) { fem.AddFixedBC(n, dof_U); fem.AddFixedBC(n, dof_V); fem.AddFixedBC(n, dof_W); }
                else if (strcmp(sz, "sxy" ) == 0) { fem.AddFixedBC(n, dof_SX); fem.AddFixedBC(n, dof_SY); }
                else if (strcmp(sz, "syz" ) == 0) { fem.AddFixedBC(n, dof_SY); fem.AddFixedBC(n, dof_SZ); }
                else if (strcmp(sz, "sxz" ) == 0) { fem.AddFixedBC(n, dof_SX); fem.AddFixedBC(n, dof_SZ); }
                else if (strcmp(sz, "sxyz") == 0) { fem.AddFixedBC(n, dof_SX); fem.AddFixedBC(n, dof_SY); fem.AddFixedBC(n, dof_SZ); }
                else if (strcmp(sz, "wxy" ) == 0) { fem.AddFixedBC(n, dof_WX); fem.AddFixedBC(n, dof_WY); }
                else if (strcmp(sz, "wyz" ) == 0) { fem.AddFixedBC(n, dof_WY); fem.AddFixedBC(n, dof_WZ); }
                else if (strcmp(sz, "wxz" ) == 0) { fem.AddFixedBC(n, dof_WX); fem.AddFixedBC(n, dof_WZ); }
                else if (strcmp(sz, "wxyz") == 0) { fem.AddFixedBC(n, dof_WX); fem.AddFixedBC(n, dof_WY); fem.AddFixedBC(n, dof_WZ); }
				else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);
			}
		}
	}
	else
	{
		// Read the fixed nodes
		++tag;
		do
		{
			int n = ReadNodeID(tag);
			const char* sz = tag.AttributeValue("bc");

			int ndof = dofs.GetDOF(sz);
			if (ndof >= 0) fem.AddFixedBC(n, ndof);
			else
			{
				// The supported fixed BC strings don't quite follow the dof naming convention.
				// For now, we'll check these BC explicitly, but I want to get rid of this in the future.
				if      (strcmp(sz, "xy"  ) == 0) { fem.AddFixedBC(n, dof_X ); fem.AddFixedBC(n, dof_Y ); }
				else if (strcmp(sz, "yz"  ) == 0) { fem.AddFixedBC(n, dof_Y ); fem.AddFixedBC(n, dof_Z ); }
				else if (strcmp(sz, "xz"  ) == 0) { fem.AddFixedBC(n, dof_X ); fem.AddFixedBC(n, dof_Z ); }
				else if (strcmp(sz, "xyz" ) == 0) { fem.AddFixedBC(n, dof_X ); fem.AddFixedBC(n, dof_Y ); fem.AddFixedBC(n, dof_Z); }
				else if (strcmp(sz, "uv"  ) == 0) { fem.AddFixedBC(n, dof_U ); fem.AddFixedBC(n, dof_V ); }
				else if (strcmp(sz, "vw"  ) == 0) { fem.AddFixedBC(n, dof_V ); fem.AddFixedBC(n, dof_W ); }
				else if (strcmp(sz, "uw"  ) == 0) { fem.AddFixedBC(n, dof_U ); fem.AddFixedBC(n, dof_W ); }
				else if (strcmp(sz, "uvw" ) == 0) { fem.AddFixedBC(n, dof_U ); fem.AddFixedBC(n, dof_V ); fem.AddFixedBC(n, dof_W); }
	            else if (strcmp(sz, "sxy" ) == 0) { fem.AddFixedBC(n, dof_SX); fem.AddFixedBC(n, dof_SY); }
		        else if (strcmp(sz, "syz" ) == 0) { fem.AddFixedBC(n, dof_SY); fem.AddFixedBC(n, dof_SZ); }
				else if (strcmp(sz, "sxz" ) == 0) { fem.AddFixedBC(n, dof_SX); fem.AddFixedBC(n, dof_SZ); }
				else if (strcmp(sz, "sxyz") == 0) { fem.AddFixedBC(n, dof_SX); fem.AddFixedBC(n, dof_SY); fem.AddFixedBC(n, dof_SZ); }
                else if (strcmp(sz, "wxy" ) == 0) { fem.AddFixedBC(n, dof_WX); fem.AddFixedBC(n, dof_WY); }
                else if (strcmp(sz, "wyz" ) == 0) { fem.AddFixedBC(n, dof_WY); fem.AddFixedBC(n, dof_WZ); }
                else if (strcmp(sz, "wxz" ) == 0) { fem.AddFixedBC(n, dof_WX); fem.AddFixedBC(n, dof_WZ); }
                else if (strcmp(sz, "wxyz") == 0) { fem.AddFixedBC(n, dof_WX); fem.AddFixedBC(n, dof_WY); fem.AddFixedBC(n, dof_WZ); }
				else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);
			}
			++tag;
		}
		while (!tag.isend());
	}
}


//-----------------------------------------------------------------------------
void FEBioBoundarySection2::ParseBCFix(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	DOFS& dofs = fem.GetDOFS();
	FEMesh& mesh = fem.GetMesh();
	int NN = mesh.Nodes();

	// get the required bc attribute
	char szbc[8];
	strcpy(szbc, tag.AttributeValue("bc"));

	// process the bc string
	vector<int> bc;

	int ndof = dofs.GetDOF(szbc);
	if (ndof >= 0) bc.push_back(ndof);
	else
	{
		// get the DOF indices
		const int dof_X = fem.GetDOFIndex("x");
		const int dof_Y = fem.GetDOFIndex("y");
		const int dof_Z = fem.GetDOFIndex("z");
		const int dof_U = fem.GetDOFIndex("u");
		const int dof_V = fem.GetDOFIndex("v");
		const int dof_W = fem.GetDOFIndex("w");
		const int dof_SX = fem.GetDOFIndex("sx");
		const int dof_SY = fem.GetDOFIndex("sy");
		const int dof_SZ = fem.GetDOFIndex("sz");
        const int dof_WX = fem.GetDOFIndex("wx");
        const int dof_WY = fem.GetDOFIndex("wy");
        const int dof_WZ = fem.GetDOFIndex("wz");

		// The supported fixed BC strings don't quite follow the dof naming convention.
		// For now, we'll check these BC explicitly, but I want to get rid of this in the future.
		if      (strcmp(szbc, "xy"  ) == 0) { bc.push_back(dof_X ); bc.push_back(dof_Y); }
		else if (strcmp(szbc, "yz"  ) == 0) { bc.push_back(dof_Y ); bc.push_back(dof_Z); }
		else if (strcmp(szbc, "xz"  ) == 0) { bc.push_back(dof_X ); bc.push_back(dof_Z); }
		else if (strcmp(szbc, "xyz" ) == 0) { bc.push_back(dof_X ); bc.push_back(dof_Y); bc.push_back(dof_Z); }
		else if (strcmp(szbc, "uv"  ) == 0) { bc.push_back(dof_U ); bc.push_back(dof_V); }
		else if (strcmp(szbc, "vw"  ) == 0) { bc.push_back(dof_V ); bc.push_back(dof_W); }
		else if (strcmp(szbc, "uw"  ) == 0) { bc.push_back(dof_U ); bc.push_back(dof_W); }
		else if (strcmp(szbc, "uvw" ) == 0) { bc.push_back(dof_U ); bc.push_back(dof_V); bc.push_back(dof_W); }
	    else if (strcmp(szbc, "sxy" ) == 0) { bc.push_back(dof_SX); bc.push_back(dof_SY); }
		else if (strcmp(szbc, "syz" ) == 0) { bc.push_back(dof_SY); bc.push_back(dof_SZ); }
		else if (strcmp(szbc, "sxz" ) == 0) { bc.push_back(dof_SX); bc.push_back(dof_SZ); }
		else if (strcmp(szbc, "sxyz") == 0) { bc.push_back(dof_SX); bc.push_back(dof_SY); bc.push_back(dof_SZ); }
        else if (strcmp(szbc, "wxy" ) == 0) { bc.push_back(dof_WX); bc.push_back(dof_WY); }
        else if (strcmp(szbc, "wyz" ) == 0) { bc.push_back(dof_WY); bc.push_back(dof_WZ); }
        else if (strcmp(szbc, "wxz" ) == 0) { bc.push_back(dof_WX); bc.push_back(dof_WZ); }
        else if (strcmp(szbc, "wxyz") == 0) { bc.push_back(dof_WX); bc.push_back(dof_WY); bc.push_back(dof_WZ); }
		else if (strcmp(szbc, "xyzuvw") == 0)
		{
			bc.push_back(dof_X); bc.push_back(dof_Y); bc.push_back(dof_Z);
			bc.push_back(dof_U); bc.push_back(dof_V); bc.push_back(dof_W);
		}
        else if (strcmp(szbc, "xyzsxyz") == 0)
        {
            bc.push_back(dof_X); bc.push_back(dof_Y); bc.push_back(dof_Z);
            bc.push_back(dof_SX); bc.push_back(dof_SY); bc.push_back(dof_SZ);
        }
		else
		{
			// see if this is a comma seperated list
			if (dofs.ParseDOFString(szbc, bc) == false)
				throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);
		}
	}

	if (bc.empty()) throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);
	int nbc = (int)bc.size();

	// create the fixed BC's
	vector<FEFixedBC*> pbc(nbc);
	for (int i=0; i<nbc; ++i)
	{
		FEFixedBC* pbci = dynamic_cast<FEFixedBC*>(fecore_new<FEBoundaryCondition>("fix", &fem));
		pbci->SetDOF(bc[i]);
		pbc[i] = pbci;

		// add it to the model
		GetBuilder()->AddFixedBC(pbci);
	}

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
		for (int j=0; j<nbc; ++j) pbc[j]->AddNodes(ns);
	}
	else
	{
		// Read the fixed nodes
		++tag;
		do
		{
			int n = ReadNodeID(tag);
			if ((n<0) || (n >= NN)) throw XMLReader::InvalidAttributeValue(tag, "id");
			for (int j=0; j<nbc; ++j) pbc[j]->AddNode(n);
			++tag;
		}
		while (!tag.isend());
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection25::ParseBCFix(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	DOFS& dofs = fem.GetDOFS();
	FEMesh& mesh = fem.GetMesh();

	// make sure the tag is a leaf
	if (tag.isempty() == false) throw XMLReader::InvalidValue(tag);

	// get the required bc attribute
	char szbc[64] = {0};
	strcpy(szbc, tag.AttributeValue("bc"));

	// process the bc string
	vector<int> bc;
	dofs.ParseDOFString(szbc, bc);

	// check the list
	if (bc.empty()) throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);
	int nbc = (int)bc.size();
	for (int i=0; i<nbc; ++i) if (bc[i] == -1) throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

	// get the nodeset
	const char* szset = tag.AttributeValue("node_set");
	map<string, FENodeSet*>::iterator nset = m_NodeSet.find(szset);
	if (nset == m_NodeSet.end()) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);
	FENodeSet* nodeSet = (*nset).second;

	// create the fixed BC's
	for (int i=0; i<nbc; ++i)
	{
		FEFixedBC* pbc = dynamic_cast<FEFixedBC*>(fecore_new<FEBoundaryCondition>("fix", &fem));
		pbc->SetDOF(bc[i]);
		pbc->AddNodes(*nodeSet);

		// add it to the model
		GetBuilder()->AddFixedBC(pbc);
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseBCPrescribe(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	DOFS& dofs = fem.GetDOFS();

	// see if this tag defines a set
	const char* szset = tag.AttributeValue("set", true);
	if (szset)
	{
		// Find the set
		FENodeSet* ps = mesh.FindNodeSet(szset);
		if (ps == 0) throw XMLReader::InvalidAttributeValue(tag, "set", szset);

		// get the bc attribute
		const char* sz = tag.AttributeValue("bc");

		// find the dof index from its symbol
		int bc = dofs.GetDOF(sz);
		if (bc == -1) 
		{
			// the temperature degree of freedom was renamed
			// for backward compatibility we need to check for it
			if (strcmp(sz, "t") == 0) bc = dofs.GetDOF("T");
			throw XMLReader::InvalidAttributeValue(tag, "bc", sz);
		}

		// get the lc attribute
		sz = tag.AttributeValue("lc");
		int lc = atoi(sz);

		// make sure this tag is a leaf
		if (tag.isleaf() == false) throw XMLReader::InvalidValue(tag);

		// get the scale factor
		double s = 1;
		value(tag, s);

		// create the bc
		FEPrescribedDOF* pdc = dynamic_cast<FEPrescribedDOF*>(fecore_new<FEBoundaryCondition>("prescribe", &fem));
		pdc->SetScale(s, lc).SetDOF(bc);

		// add this boundary condition to the current step
		GetBuilder()->AddPrescribedBC(pdc);

		// add nodes in the nodeset
		pdc->AddNodes(*ps);
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
			int n = ReadNodeID(tag);
			const char* sz = tag.AttributeValue("bc");

			// get the dof index from its symbol
			int bc = dofs.GetDOF(sz);
			if (bc == -1) throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

			sz = tag.AttributeValue("lc");
			int lc = atoi(sz)-1;

			double scale;
			tag.value(scale);

			FEPrescribedDOF* pdc = dynamic_cast<FEPrescribedDOF*>(fecore_new<FEBoundaryCondition>("prescribe", &fem));
			pdc->SetDOF(bc).SetScale(scale, lc).SetRelativeFlag(br);
			pdc->AddNode(n);

			// add this boundary condition to the current step
			GetBuilder()->AddPrescribedBC(pdc);

			// next tag
			++tag;
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection2::ParseBCPrescribe(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	DOFS& dofs = fem.GetDOFS();
	int NN = mesh.Nodes();

	// count how many prescibed nodes there are
	int ndis = tag.children();

	// determine whether prescribed BC is relative or absolute
	bool br = false;
	const char* sztype = tag.AttributeValue("type",true);
	if (sztype && strcmp(sztype, "relative") == 0) br = true;

	// get the BC
	const char* sz = tag.AttributeValue("bc");
	int bc = dofs.GetDOF(sz);
	if (bc == -1) 
	{
		// the temperature degree of freedom was renamed
		// for backward compatibility we need to check for it
		if (strcmp(sz, "t") == 0) bc = dofs.GetDOF("T");
		else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);
	}

	// get the load curve number
	int lc = -1;
	sz = tag.AttributeValue("lc", true);
	if (sz) { lc = atoi(sz) - 1; }

	// see if the scale attribute is defined
	double scale = 1.0;
	tag.AttributeValue("scale", scale, true);

	// if the value is a leaf, the scale factor is defined as the value
	if (tag.isleaf())
	{
		tag.value(scale);
	}

	// create a prescribed bc
	FEPrescribedDOF* pdc = dynamic_cast<FEPrescribedDOF*>(fecore_new<FEBoundaryCondition>("prescribe", &fem));
	pdc->SetDOF(bc).SetScale(scale, lc).SetRelativeFlag(br);

	// add this boundary condition to the current step
	GetBuilder()->AddPrescribedBC(pdc);

	// see if there is a set defined
	const char* szset = tag.AttributeValue("set", true);
	if (szset)
	{
		// make sure this is a leaf tag
		if (tag.isleaf() == false) throw XMLReader::InvalidValue(tag);

		// find the node set
		FENodeSet* pns = mesh.FindNodeSet(szset);
		if (pns == 0) throw XMLReader::InvalidAttributeValue(tag, "set", szset);

		// add the nodes
		pdc->AddNodes(*pns);
	}
	else
	{
		// read the prescribed data
		++tag;
		for (int i=0; i<ndis; ++i)
		{
			// get the node ID
			int n = ReadNodeID(tag);
			value(tag, scale);

			pdc->AddNode(n, scale);
			++tag;
		}
	}
}

//-----------------------------------------------------------------------------
//! In version 2.5 all prescribed boundary conditions use a node set to define
//! the list of nodes to which the bc is applied.
void FEBioBoundarySection25::ParseBCPrescribe(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	DOFS& dofs = fem.GetDOFS();
	int NN = mesh.Nodes();

	// get the BC
	const char* sz = tag.AttributeValue("bc");
	int bc = dofs.GetDOF(sz);
	if (bc == -1) throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

	// get the node set (if defined)
	const char* szset = tag.AttributeValue("node_set");
	map<string,FENodeSet*>::iterator nset = m_NodeSet.find(szset);
	if (nset == m_NodeSet.end()) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);
	FENodeSet* nodeSet = (*nset).second;

	// create a prescribed bc
	FEPrescribedDOF* pdc = dynamic_cast<FEPrescribedDOF*>(fecore_new<FEBoundaryCondition>("prescribe", &fem));
	pdc->SetDOF(bc);
	pdc->AddNodes(*nodeSet);

	// add this boundary condition to the current step
	GetBuilder()->AddPrescribedBC(pdc);

	// Read the parameter list
	FEParameterList& pl = pdc->GetParameterList();
	ReadParameterList(tag, pl);
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection25::ParseBC(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the type string
	const char* sztype = tag.AttributeValue("type");

	// create the boundary condition
	FEPrescribedBC* pdc = dynamic_cast<FEPrescribedBC*>(fecore_new<FEBoundaryCondition>(sztype, &fem));
	if (pdc == 0) throw XMLReader::InvalidTag(tag);

	// get the node set
	const char* szset = tag.AttributeValue("node_set", true);
	if (szset)
	{
		FENodeSet* nodeSet = mesh.FindNodeSet(szset);
		if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

		// add the nodes to the BC
		pdc->AddNodes(*nodeSet);

		// Read the parameter list
		FEParameterList& pl = pdc->GetParameterList();
		ReadParameterList(tag, pl);
	}
	else
	{
		// if a node set is not defined, see if a surface is defined
		szset = tag.AttributeValue("surface");
		FEFacetSet* set = mesh.FindFacetSet(szset);
		FESurface* surf = new FESurface(set);

		// Read the parameter list (before setting the surface)
		FEParameterList& pl = pdc->GetParameterList();
		ReadParameterList(tag, pl);

		// add the surface nodes
		pdc->AddNodes(*surf);

		// don't forget to cleanup
		delete surf;
	}

	// add this boundary condition to the current step
	GetBuilder()->AddPrescribedBC(pdc);
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseSpringSection(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// determine the spring type
	const char* szt = tag.AttributeValue("type", true);
	if (szt == 0) szt = "linear";
	FEDiscreteMaterial* pm = dynamic_cast<FEDiscreteMaterial*>(fecore_new<FEMaterial>(szt, &fem));
	if (pm == 0) throw XMLReader::InvalidAttributeValue(tag, "type", szt);

	// create a new spring "domain"
	FECoreKernel& febio = FECoreKernel::GetInstance();
	FE_Element_Spec spec;
	spec.eclass = FE_ELEM_TRUSS;
	spec.eshape = ET_TRUSS2;
	spec.etype  = FE_DISCRETE;
	FEDiscreteDomain* pd = dynamic_cast<FEDiscreteDomain*>(febio.CreateDomain(spec, &mesh, pm));
	mesh.AddDomain(pd);

	pd->Create(1, FE_DISCRETE);
	FEDiscreteElement& de = pd->Element(0);
	de.SetID(++GetBuilder()->m_maxid);
	
	// add a new material for each spring
	fem.AddMaterial(pm);
	pm->SetID(fem.Materials());
	pd->SetMatID(fem.Materials()-1);

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
			if (ReadParameter(tag, pl) == 0)
			{
				throw XMLReader::InvalidTag(tag);
			}
		}
		++tag;
	}
	while (!tag.isend());

	pd->CreateMaterialPointData();
}

//-----------------------------------------------------------------------------
//! Parse the linear constraints section of the xml input file
//! This section is a subsection of the Boundary section

void FEBioBoundarySection::ParseConstraints(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	DOFS& dofs = fem.GetDOFS();

	// make sure there is a constraint defined
	if (tag.isleaf()) return;

	FEModelBuilder* feb = GetBuilder();

	// read the master node
	int nodeID;
	tag.AttributeValue("node", nodeID);
	int masterNode = feb->FindNodeFromID(nodeID);

	// get the dofs
	const char* szbc = tag.AttributeValue("bc");
	vector<int> dofList;
	dofs.ParseDOFString(szbc, dofList);
	int ndofs = (int) dofList.size();

	// allocate linear constraints
	vector<FELinearConstraint> LC(ndofs, FELinearConstraint(&fem));
	for (int i=0; i<ndofs; ++i)
	{
		int dof = dofList[i];
		if (dof < 0) throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);
		LC[i].master.dof = dof;
		LC[i].master.node = masterNode;
	}

	// read the slave nodes
	++tag;
	do
	{
		FELinearConstraint::DOF dof;
		if (tag == "node")
		{
			// get the node
			int slaveNode = ReadNodeID(tag);

			// get the dof
			// (if ommitted we take the master dof)
			int slaveDOF = -1;
			const char* szbc = tag.AttributeValue("bc", true);
			if (szbc)
			{
				slaveDOF = dofs.GetDOF(szbc);
				if (slaveDOF < 0) throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);
			}

			// get the coefficient
			double val;
			tag.value(val);

			// add it to the list
			for (int i=0; i<ndofs; ++i)
			{
				dof.node = slaveNode;
				dof.dof  = (slaveDOF < 0 ? LC[i].master.dof : slaveDOF);
				dof.val  = val;
				LC[i].slave.push_back(dof);
			}
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	// add the linear constraint to the system
	for (int i=0; i<ndofs; ++i)
	{
		fem.GetLinearConstraintManager().AddLinearConstraint(LC[i]);

		GetBuilder()->AddComponent(&LC[i]);
	}
}

void FEBioBoundarySection25::ParseMergeConstraint(XMLTag& tag)
{
	// make sure this is an empty tag
	if (tag.isempty() == false) throw XMLReader::InvalidValue(tag);

	// get the dofs
	const char* szbc = tag.AttributeValue("bc");

	// get the surface pair name
	const char* szsp = tag.AttributeValue("surface_pair");

	// get the dof list
	FEModel& fem = *GetFEModel();
	DOFS& dof = fem.GetDOFS();
	vector<int> dofs;
	dof.ParseDOFString(szbc, dofs);
	if (dofs.empty()) throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

	// get the surfaces
	FEMesh& mesh = fem.GetMesh();
	FESurfacePair* sp = mesh.FindSurfacePair(szsp);
	if (sp == 0) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szsp);

	// merge the interfaces
	FEMergedConstraint merge(fem);
	if (merge.Merge(sp->GetMasterSurface(), sp->GetSlaveSurface(), dofs) == false)
		throw XMLReader::InvalidTag(tag);
}

void FEBioBoundarySection25::ParsePeriodicLinearConstraint(XMLTag& tag)
{
	FEModel* fem = GetFEModel();
	FEMesh& mesh = fem->GetMesh();
	FEPeriodicLinearConstraint plc;

	FEModelBuilder* feb = GetBuilder();

	++tag;
	do
	{
		if (tag == "constrain")
		{
			const char* sz = tag.AttributeValue("surface_pair", true);
			if (sz)
			{
				FESurfacePair* spair = mesh.FindSurfacePair(sz);
				if (spair == 0) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", sz);

				FESurface* ms = new FESurface(spair->GetMasterSurface()); feb->BuildSurface(*ms, *spair->GetMasterSurface());
				FESurface* ss = new FESurface(spair->GetSlaveSurface ()); feb->BuildSurface(*ss, *spair->GetSlaveSurface());
				plc.AddNodeSetPair(ms->GetNodeSet(), ss->GetNodeSet());
			}
			else throw XMLReader::MissingAttribute(tag, "surface_pair");
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	// generate the linear constraints
	plc.GenerateConstraints(fem);

	// don't forget to activate

}

void FEBioBoundarySection25::ParsePeriodicLinearConstraint2O(XMLTag& tag)
{
	FEModel* fem = GetFEModel();
	FEMesh& mesh = fem->GetMesh();
	FEPeriodicLinearConstraint2O plc;

	FEModelBuilder* feb = GetBuilder();

	++tag;
	do
	{
		if (tag == "constrain")
		{
			const char* sz = tag.AttributeValue("surface_pair");
			if (sz)
			{
				FESurfacePair* spair = mesh.FindSurfacePair(sz);
				if (spair == 0) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", sz);

				FESurface* ms = new FESurface(spair->GetMasterSurface()); feb->BuildSurface(*ms, *spair->GetMasterSurface());
				FESurface* ss = new FESurface(spair->GetSlaveSurface()); feb->BuildSurface(*ss, *spair->GetSlaveSurface());
				plc.AddNodeSetPair(ms->GetNodeSet(), ss->GetNodeSet());
			}
			else throw XMLReader::MissingAttribute(tag, "surface_pair");
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	// generate the linear constraints
	if (plc.GenerateConstraints(fem) == false)
	{
		throw XMLReader::InvalidTag(tag);
	}
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection::ParseContactInterface(XMLTag& tag, FESurfacePairConstraint* pci)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	// get the parameter list
	FEParameterList& pl = pci->GetParameterList();

	// read the parameters
	++tag;
	do
	{
		if (ReadParameter(tag, pl) == false)
		{
			if (tag == "surface")
			{
				const char* sztype = tag.AttributeValue("type");
				int ntype = 0;
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
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidSystem& rigid = *fem.GetRigidSystem();
	FEMesh& m = fem.GetMesh();

	FEModelBuilder* feb = GetBuilder();

	// get the type attribute
	const char* szt = tag.AttributeValue("type");

	// Not all contact interfaces can be parsed automatically.
	// First, check all these special cases.
	if (strcmp(szt, "rigid_wall") == 0)
	{
		// --- R I G I D   W A L L   I N T E R F A C E ---

		FERigidWallInterface* ps = dynamic_cast<FERigidWallInterface*>(fecore_new<FESurfacePairConstraint>(szt, GetFEModel()));
		if (ps)
		{
			fem.AddSurfacePairConstraint(ps);

			++tag;
			do
			{
				if (tag == "plane")
				{
					// In old formats, the load curve was set in the "plane" property.
					// Now, we need to map this load curve to the "offset" parameter.
					const char* sz = tag.AttributeValue("lc", true);
					if (sz)
					{
						int nlc = atoi(sz)-1;
						FEParameterList& pl = ps->GetParameterList();
						FEParam& p = *pl.FindFromName("offset");
						p.SetLoadCurve(nlc, 1.0);
					}
				}

				if (ReadParameter(tag, ps) == false)
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
		int id, rb, rbp = -1;
		FERigidNodeSet* prn = 0;
		for (int i=0; i<nrn; ++i)
		{
			id = atoi(tag.AttributeValue("id"))-1;
			rb = atoi(tag.AttributeValue("rb"))-1;

			if ((prn == 0) || (rb != rbp))
			{
				prn = new FERigidNodeSet(&fem);
				prn->SetRigidID(rb);

				// the default shell bc depends on the shell formulation
				prn->SetShellBC(feb->m_default_shell == OLD_SHELL ? FERigidNodeSet::HINGED_SHELL : FERigidNodeSet::CLAMPED_SHELL);

				rigid.AddRigidNodeSet(prn);

				feb->AddComponent(prn);
				rbp = rb;
			}
			prn->AddNode(id);

			++tag;
		}
	}
	else if (strcmp(szt, "linear constraint") == 0)
	{
		FEModel& fem = *GetFEModel();

		// make sure there is a constraint defined
		if (tag.isleaf()) return;

		// create a new linear constraint manager
		FELinearConstraintSet* pLCS = dynamic_cast<FELinearConstraintSet*>(fecore_new<FENLConstraint>(szt, GetFEModel()));
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
						dof.node = ReadNodeID(tag);

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
			else if (ReadParameter(tag, pLCS) == false)
			{
				throw XMLReader::InvalidTag(tag);
			}
			++tag;
		}
		while (!tag.isend());
	}
	else
	{
		// If we get here, we try to create a contact interface
		// using the FEBio kernel. 
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fecore_new<FESurfacePairConstraint>(szt, GetFEModel()));
		if (pci)
		{
			// add it to the model
			GetBuilder()->AddContactInterface(pci);

			// parse the interface
			ParseContactInterface(tag, pci);
		}
		else 
		{
			// some nonlinear constraints are also defined in the Contact section, so let's try that next.
			// TODO: These are mostly rigid constraints. Therefore, I would like to move this elsewhere (maybe in the new Rigid section?)
			FENLConstraint* pnlc = fecore_new<FENLConstraint>(szt, GetFEModel());
			if (pnlc)
			{
				ReadParameterList(tag, pnlc);
				fem.AddNonlinearConstraint(pnlc);
			}
			else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
		}
	}
}

//-----------------------------------------------------------------------------
// Rigid node sets are defined in the Boundary section since version 2.5
// (Used to be defined in the Contact section)
void FEBioBoundarySection25::ParseBCRigid(XMLTag& tag)
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FEMesh& mesh = fem.GetMesh();
	FERigidSystem& rigid = *fem.GetRigidSystem();
	FEModelBuilder* feb = GetBuilder();
	int NMAT = fem.Materials();

	// get the rigid body material ID
	int rb = -1;
	tag.AttributeValue("rb", rb);
	rb -= 1;

	// make sure we have a valid rigid body reference
	if ((rb < 0)||(rb>=NMAT)) throw XMLReader::InvalidAttributeValue(tag, "rb", tag.AttributeValue("rb"));

	// get the nodeset
	const char* szset = tag.AttributeValue("node_set");
	FENodeSet* nodeSet = mesh.FindNodeSet(szset);
	if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

	// create new rigid node set
	FERigidNodeSet* prn = new FERigidNodeSet(&fem);

	// the default shell bc depends on the shell formulation
	prn->SetShellBC(feb->m_default_shell == OLD_SHELL ? FERigidNodeSet::HINGED_SHELL : FERigidNodeSet::CLAMPED_SHELL);

	prn->SetRigidID(rb);
	prn->SetNodeSet(*nodeSet);

	rigid.AddRigidNodeSet(prn);
	
	// add it to the current step
	GetBuilder()->AddComponent(prn);

	// read the parameter list
	ReadParameterList(tag, prn);
}

//-----------------------------------------------------------------------------
// The rigid body "constraints" are moved to the Boundary section in 2.5
void FEBioBoundarySection25::ParseRigidBody(XMLTag& tag)
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidSystem& rigid = *fem.GetRigidSystem();

	const char* szm = tag.AttributeValue("mat");
	assert(szm);

	// get the material ID
	int nmat = atoi(szm);
	if ((nmat <= 0) || (nmat > fem.Materials())) throw XMLReader::InvalidAttributeValue(tag, "mat", szm);

	// make sure this is a valid rigid material
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(nmat-1));
	if (pm == 0) throw XMLReader::InvalidAttributeValue(tag, "mat", szm);

	++tag;
	do
	{
		if (tag == "prescribed")
		{
			// get the dof
			int bc = -1;
			const char* szbc = tag.AttributeValue("bc");
			if      (strcmp(szbc, "x") == 0) bc = 0;
			else if (strcmp(szbc, "y") == 0) bc = 1;
			else if (strcmp(szbc, "z") == 0) bc = 2;
			else if (strcmp(szbc, "Rx") == 0) bc = 3;
			else if (strcmp(szbc, "Ry") == 0) bc = 4;
			else if (strcmp(szbc, "Rz") == 0) bc = 5;
			else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

			// get the loadcurve
			const char* szlc = tag.AttributeValue("lc");
			int lc = atoi(szlc) - 1;

			// get the (optional) type attribute
			bool brel = false;
			const char* szrel = tag.AttributeValue("type", true);
			if (szrel)
			{
				if      (strcmp(szrel, "relative" ) == 0) brel = true;
				else if (strcmp(szrel, "absolute" ) == 0) brel = false;
				else throw XMLReader::InvalidAttributeValue(tag, "type", szrel);
			}

			// create the rigid displacement constraint
			FERigidBodyDisplacement* pDC = static_cast<FERigidBodyDisplacement*>(fecore_new<FEBoundaryCondition>("rigid_prescribed", &fem));
			pDC->id = nmat;
			pDC->bc = bc;
			pDC->lc = lc;
			pDC->brel = brel;
			value(tag, pDC->sf);
			rigid.AddPrescribedBC(pDC);

			// add this boundary condition to the current step
			GetBuilder()->AddComponent(pDC);
		}
		else if (tag == "force")
		{
			// get the dof
			int bc = -1;
			const char* szbc = tag.AttributeValue("bc");
			if      (strcmp(szbc, "x") == 0) bc = 0;
			else if (strcmp(szbc, "y") == 0) bc = 1;
			else if (strcmp(szbc, "z") == 0) bc = 2;
			else if (strcmp(szbc, "Rx") == 0) bc = 3;
			else if (strcmp(szbc, "Ry") == 0) bc = 4;
			else if (strcmp(szbc, "Rz") == 0) bc = 5;
			else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

			// get the type
			int ntype = 0;
			bool bfollow = false;
			const char* sztype = tag.AttributeValue("type", true);
			if (sztype)
			{
				if (strcmp(sztype, "ramp") == 0) ntype = 1;
				else if (strcmp(sztype, "follow") == 0) bfollow = true;
				else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
			}

			// get the loadcurve
			const char* szlc = tag.AttributeValue("lc", true);
			int lc = -1;
			if (szlc) lc = atoi(szlc) - 1;

			// make sure there is a loadcurve for type=0 forces
			if ((ntype == 0)&&(lc==-1)) throw XMLReader::MissingAttribute(tag, "lc");

			// create the rigid body force
			FERigidBodyForce* pFC = static_cast<FERigidBodyForce*>(fecore_new<FEModelLoad>(FEBC_ID, "rigid_force",  &fem));
			pFC->m_ntype = ntype;
			pFC->id = nmat;
			pFC->bc = bc;
			pFC->lc = lc;
			pFC->m_bfollow = bfollow;
			value(tag, pFC->sf);

			// add it to the model
			GetBuilder()->AddModelLoad(pFC);
		}
		else if (tag == "fixed")
		{
			// get the dof
			int bc = -1;
			const char* szbc = tag.AttributeValue("bc");
			if      (strcmp(szbc, "x") == 0) bc = 0;
			else if (strcmp(szbc, "y") == 0) bc = 1;
			else if (strcmp(szbc, "z") == 0) bc = 2;
			else if (strcmp(szbc, "Rx") == 0) bc = 3;
			else if (strcmp(szbc, "Ry") == 0) bc = 4;
			else if (strcmp(szbc, "Rz") == 0) bc = 5;
			else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

			// create the fixed dof
			FERigidBodyFixedBC* pBC = static_cast<FERigidBodyFixedBC*>(fecore_new<FEBoundaryCondition>(FEBC_ID, "rigid_fixed",  &fem));
			pBC->id = nmat;
			pBC->bc = bc;
			rigid.AddFixedBC(pBC);

			// add this boundary condition to the current step
			GetBuilder()->AddComponent(pBC);
		}
		else if (tag == "initial_velocity")
		{
			// get the initial velocity
			vec3d v;
			value(tag, v);

			// create the initial condition
			FERigidBodyVelocity* pic = new FERigidBodyVelocity(&fem);
			pic->m_rid = nmat;
			pic->m_vel = v;
			rigid.AddInitialVelocity(pic);

			// add this initial condition to the current step
			GetBuilder()->AddComponent(pic);
		}
		else if (tag == "initial_angular_velocity")
		{
			// get the initial angular velocity
			vec3d w;
			value(tag, w);

			// create the initial condition
			FERigidBodyAngularVelocity* pic = new FERigidBodyAngularVelocity(&fem);
			pic->m_rid = nmat;
			pic->m_w = w;
			rigid.AddInitialAngularVelocity(pic);

			// add this initial condition to the current step
			GetBuilder()->AddComponent(pic);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
