/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FEBioBoundarySection.h"
#include <FECore/FEModel.h>
#include <FECore/FEDiscreteMaterial.h>
#include <FECore/FEDiscreteDomain.h>
#include <FECore/FEAugLagLinearConstraint.h>
#include <FECore/FEPrescribedDOF.h>
#include <FECore/FEFixedBC.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FEPeriodicLinearConstraint.h>
//#include <FEBioRVE/FEPeriodicLinearConstraint2O.h>
#include <FECore/FEMergedConstraint.h>
#include <FEBioMech/FEMechModel.h>
#include <FEBioMech/FERigidMaterial.h>
#include <FECore/FEFacetSet.h>
#include <FECore/log.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FEInitialCondition.h>

//---------------------------------------------------------------------------------
void FEBioBoundarySection::BuildNodeSetMap()
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	m_NodeSet.clear();
	for (int i = 0; i<m.NodeSets(); ++i)
	{
		FENodeSet* nsi = m.NodeSet(i);
		m_NodeSet[nsi->GetName()] = nsi;
	}
}

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
				int nn = pe->GetFace(nf[1]-1, ne);
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
void FEBioBoundarySection::AddFixedBC(FENodeSet* set, int dof)
{
	FEModel* fem = GetFEModel();
	FEFixedBC* bc = new FEFixedBC(fem);
	bc->SetDOFList(dof);
	bc->SetNodeSet(set);
	fem->AddBoundaryCondition(bc);
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
    const int dof_GX = fem.GetDOFIndex("gx");
    const int dof_GY = fem.GetDOFIndex("gy");
    const int dof_GZ = fem.GetDOFIndex("gz");

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
		int ndof = dofs.GetDOF(sz);
		if (ndof >= 0) AddFixedBC(ps, ndof);
		else
		{
			// The supported fixed BC strings don't quite follow the dof naming convention.
			// For now, we'll check these BC explicitly, but I want to get rid of this in the future.
			if      (strcmp(sz, "xy"  ) == 0) { AddFixedBC(ps, dof_X ); AddFixedBC(ps, dof_Y); }
			else if (strcmp(sz, "yz"  ) == 0) { AddFixedBC(ps, dof_Y ); AddFixedBC(ps, dof_Z); }
			else if (strcmp(sz, "xz"  ) == 0) { AddFixedBC(ps, dof_X ); AddFixedBC(ps, dof_Z); }
			else if (strcmp(sz, "xyz" ) == 0) { AddFixedBC(ps, dof_X ); AddFixedBC(ps, dof_Y); AddFixedBC(ps, dof_Z); }
			else if (strcmp(sz, "uv"  ) == 0) { AddFixedBC(ps, dof_U ); AddFixedBC(ps, dof_V); }
			else if (strcmp(sz, "vw"  ) == 0) { AddFixedBC(ps, dof_V ); AddFixedBC(ps, dof_W); }
			else if (strcmp(sz, "uw"  ) == 0) { AddFixedBC(ps, dof_U ); AddFixedBC(ps, dof_W); }
			else if (strcmp(sz, "uvw" ) == 0) { AddFixedBC(ps, dof_U ); AddFixedBC(ps, dof_V); AddFixedBC(ps, dof_W); }
            else if (strcmp(sz, "sxy" ) == 0) { AddFixedBC(ps, dof_SX); AddFixedBC(ps, dof_SY); }
            else if (strcmp(sz, "syz" ) == 0) { AddFixedBC(ps, dof_SY); AddFixedBC(ps, dof_SZ); }
            else if (strcmp(sz, "sxz" ) == 0) { AddFixedBC(ps, dof_SX); AddFixedBC(ps, dof_SZ); }
            else if (strcmp(sz, "sxyz") == 0) { AddFixedBC(ps, dof_SX); AddFixedBC(ps, dof_SY); AddFixedBC(ps, dof_SZ); }
            else if (strcmp(sz, "wxy" ) == 0) { AddFixedBC(ps, dof_WX); AddFixedBC(ps, dof_WY); }
            else if (strcmp(sz, "wyz" ) == 0) { AddFixedBC(ps, dof_WY); AddFixedBC(ps, dof_WZ); }
            else if (strcmp(sz, "wxz" ) == 0) { AddFixedBC(ps, dof_WX); AddFixedBC(ps, dof_WZ); }
            else if (strcmp(sz, "wxyz") == 0) { AddFixedBC(ps, dof_WX); AddFixedBC(ps, dof_WY); AddFixedBC(ps, dof_WZ); }
            else if (strcmp(sz, "gxy" ) == 0) { AddFixedBC(ps, dof_GX); AddFixedBC(ps, dof_GY); }
            else if (strcmp(sz, "gyz" ) == 0) { AddFixedBC(ps, dof_GY); AddFixedBC(ps, dof_GZ); }
            else if (strcmp(sz, "gxz" ) == 0) { AddFixedBC(ps, dof_GX); AddFixedBC(ps, dof_GZ); }
            else if (strcmp(sz, "gxyz") == 0) { AddFixedBC(ps, dof_GX); AddFixedBC(ps, dof_GY); AddFixedBC(ps, dof_GZ); }
			else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);
		}
	}
	else
	{
		// The format where the bc can be defined on each line is no longer supported.
		throw XMLReader::MissingAttribute(tag, "set");
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
	char szbc[32] = { 0 };
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
        const int dof_GX = fem.GetDOFIndex("gx");
        const int dof_GY = fem.GetDOFIndex("gy");
        const int dof_GZ = fem.GetDOFIndex("gz");

		// The supported fixed BC strings don't quite follow the dof naming convention.
		// For now, we'll check these BC explicitly, but I want to get rid of this in the future.
		if      (strcmp(szbc, "x"   ) == 0) { bc.push_back(dof_X ); }
		else if (strcmp(szbc, "y"   ) == 0) { bc.push_back(dof_Y ); }
		else if (strcmp(szbc, "z"   ) == 0) { bc.push_back(dof_Z ); }
		else if (strcmp(szbc, "xy"  ) == 0) { bc.push_back(dof_X ); bc.push_back(dof_Y); }
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
        else if (strcmp(szbc, "gxy" ) == 0) { bc.push_back(dof_GX); bc.push_back(dof_GY); }
        else if (strcmp(szbc, "gyz" ) == 0) { bc.push_back(dof_GY); bc.push_back(dof_GZ); }
        else if (strcmp(szbc, "gxz" ) == 0) { bc.push_back(dof_GX); bc.push_back(dof_GZ); }
        else if (strcmp(szbc, "gxyz") == 0) { bc.push_back(dof_GX); bc.push_back(dof_GY); bc.push_back(dof_GZ); }
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

	// create the fixed BC's
	FEFixedDOF* pbc = dynamic_cast<FEFixedDOF*>(fecore_new<FEBoundaryCondition>("fix", &fem));
	pbc->SetDOFS(bc);

	// add it to the model
	GetBuilder()->AddBC(pbc);

	// see if the set attribute is defined
	const char* szset = tag.AttributeValue("set", true);
	if (szset)
	{
		// make sure the tag is a leaf
		if (tag.isleaf() == false) throw XMLReader::InvalidValue(tag);

		// process the node set
		FENodeSet* pns = mesh.FindNodeSet(szset);
		if (pns == 0) throw XMLReader::InvalidAttributeValue(tag, "set", szset);

		pbc->SetNodeSet(pns);
	}
	else
	{
		FENodeSet* nset = new FENodeSet(&fem);
		fem.GetMesh().AddNodeSet(nset);

		// Read the fixed nodes
		++tag;
		do
		{
			int n = ReadNodeID(tag);
			if ((n<0) || (n >= NN)) throw XMLReader::InvalidAttributeValue(tag, "id");
			nset->Add(n);
			++tag;
		}
		while (!tag.isend());

		pbc->SetNodeSet(nset);
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

	// create the fixed BC
	FEFixedDOF* pbc = dynamic_cast<FEFixedDOF*>(fecore_new<FEBoundaryCondition>("fix", &fem));
	pbc->SetDOFS(bc);
	pbc->SetNodeSet(nodeSet);

	// add it to the model
	GetBuilder()->AddBC(pbc);
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
		pdc->SetScale(s).SetDOF(bc);

		if (lc >= 0)
		{
			FEParam* p = pdc->GetParameter("scale");
			if (p == nullptr) throw XMLReader::InvalidTag(tag);
			fem.AttachLoadController(p, lc);
		}

		// add this boundary condition to the current step
		GetBuilder()->AddBC(pdc);

		// add nodes in the nodeset
		pdc->SetNodeSet(ps);
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
			pdc->SetDOF(bc);
			pdc->SetScale(scale);
			pdc->SetRelativeFlag(br);

			FENodeSet* ps = new FENodeSet(&fem);
			ps->Add(n);
			pdc->SetNodeSet(ps);

			if (lc >= 0)
			{
				FEParam* p = pdc->GetParameter("scale");
				if (p == nullptr) throw XMLReader::InvalidTag(tag);
				fem.AttachLoadController(p, lc);
			}

			// add this boundary condition to the current step
			GetBuilder()->AddBC(pdc);

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
	pdc->SetDOF(bc);
	pdc->SetScale(scale);
	pdc->SetRelativeFlag(br);

	if (lc >= 0)
	{
		FEParam* p = pdc->GetParameter("scale");
		if (p == nullptr) throw XMLReader::InvalidTag(tag);
		fem.AttachLoadController(p, lc);
	}

	// add this boundary condition to the current step
	GetBuilder()->AddBC(pdc);

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
		pdc->SetNodeSet(pns);
	}
	else
	{
		FENodeSet* nset = new FENodeSet(&fem);
		fem.GetMesh().AddNodeSet(nset);

		// read the prescribed data
		++tag;
		for (int i=0; i<ndis; ++i)
		{
			// get the node ID
			int n = ReadNodeID(tag);
			value(tag, scale);

			// TODO: I need to create a data map for this BC and assign the values
			// to that data map

			nset->Add(n);
			++tag;
		}

		pdc->SetNodeSet(nset);
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

	// Boundary conditions can be applied to node sets or surfaces
	// depending on whether the node_set or surface attribute is defined
	FENodeSet* nodeSet = nullptr;
	FEFacetSet* facetSet = nullptr;

	// get the node set or surface
	const char* szset = tag.AttributeValue("node_set");
	if (szset) {
		map<string, FENodeSet*>::iterator nset = m_NodeSet.find(szset);
		if (nset == m_NodeSet.end()) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);
		nodeSet = (*nset).second;
	}
	else
	{
		const char* szset = tag.AttributeValue("surface");
		if (szset) {
			facetSet = mesh.FindFacetSet(szset);
			if (facetSet == nullptr) throw XMLReader::InvalidAttributeValue(tag, "surface", szset);
		}
		else throw XMLReader::MissingAttribute(tag, "node_set");
	}

	// create a prescribed bc
	FEPrescribedDOF* pdc = dynamic_cast<FEPrescribedDOF*>(fecore_new<FEBoundaryCondition>("prescribe", &fem));
	pdc->SetDOF(bc);

	// Apply either the node set or the surface
	if (nodeSet) pdc->SetNodeSet(nodeSet);
	else if (facetSet)
	{
		FENodeSet* set = new FENodeSet(&fem);
		set->Add(facetSet->GetNodeList());
		pdc->SetNodeSet(set);
	}

	// add this boundary condition to the current step
	GetBuilder()->AddBC(pdc);

	// Read the parameter list
	FEParameterList& pl = pdc->GetParameterList();

	++tag;
	do
	{
		if (ReadParameter(tag, pl, 0, 0) == false)
		{
			if (tag == "value")
			{
				feLogWarningEx((&fem), "The value parameter of the prescribed bc is deprecated.");

				// NOTE: This will only work if the scale was set to 1!!
				const char* sznodedata = tag.AttributeValue("node_data", true);
				if (sznodedata)
				{
					FEParam* pp = pdc->GetParameter("scale"); assert(pp);
					GetBuilder()->AddMappedParameter(pp, pdc, sznodedata);
				}
				else
				{
					double v;
					tag.value(v);
					pdc->SetScale(v);
				}
			}
			else throw XMLReader::InvalidTag(tag);
		}
		++tag;
	} while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection25::ParseBC(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the type string
	const char* sztype = tag.AttributeValue("type");

	// create the boundary condition
	FEBoundaryCondition* pdc = fecore_new<FEBoundaryCondition>(sztype, &fem);
	if (pdc == 0) throw XMLReader::InvalidTag(tag);

	// get the selection
	const char* szset = tag.AttributeValue("node_set", true);
	if (dynamic_cast<FENodalBC*>(pdc))
	{
		FENodalBC* pnbc = dynamic_cast<FENodalBC*>(pdc);
		FENodeSet* nodeSet = mesh.FindNodeSet(szset);
		if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

		// Set the node set
		pnbc->SetNodeSet(nodeSet);

		// Read the parameter list
		FEParameterList& pl = pdc->GetParameterList();
		ReadParameterList(tag, pl);
	}
	else if (dynamic_cast<FESurfaceBC*>(pdc))
	{
		FESurfaceBC* sbc = dynamic_cast<FESurfaceBC*>(pdc);

		// if a node set is not defined, see if a surface is defined
		szset = tag.AttributeValue("surface");
		FEFacetSet* set = mesh.FindFacetSet(szset);

		// Read the parameter list (before setting the surface)
		FEParameterList& pl = pdc->GetParameterList();
		ReadParameterList(tag, pl);

		// add the surface
		FESurface* surf = fecore_alloc(FESurface, GetFEModel());
		surf->Create(*set);
		mesh.AddSurface(surf);

		sbc->SetSurface(surf);
	}

	// add this boundary condition to the current step
	GetBuilder()->AddBC(pdc);
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

	pd->Create(1, spec);
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

	// read the parent node
	int nodeID;
	tag.AttributeValue("node", nodeID);
	int parentNode = feb->FindNodeFromID(nodeID);

	// get the dofs
	const char* szbc = tag.AttributeValue("bc");
	vector<int> dofList;
	dofs.ParseDOFString(szbc, dofList);
	int ndofs = (int) dofList.size();

	// allocate linear constraints
	vector<FELinearConstraint*> LC;
	for (int i=0; i<ndofs; ++i)
	{
		int dof = dofList[i];
		if (dof < 0) throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

		LC[i] = fecore_alloc(FELinearConstraint, &fem);
		LC[i]->SetParentDof(dof, parentNode);
	}

	// read the child nodes
	++tag;
	do
	{
		if (tag == "node")
		{
			// get the node
			int childNode = ReadNodeID(tag);

			// get the dof
			// (if ommitted we take the parent dof)
			int childDOF = -1;
			const char* szbc = tag.AttributeValue("bc", true);
			if (szbc)
			{
				childDOF = dofs.GetDOF(szbc);
				if (childDOF < 0) throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);
			}

			// get the coefficient
			double val;
			tag.value(val);

			// add it to the list
			for (int i=0; i<ndofs; ++i)
			{
				int ndof = (childDOF < 0 ? LC[i]->GetParentDof() : childDOF);
				LC[i]->AddChildDof(ndof, childNode, val);
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
		GetBuilder()->AddComponent(LC[i]);
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
	if (merge.Merge(sp->GetSecondarySurface(), sp->GetPrimarySurface(), dofs) == false)
		throw XMLReader::InvalidTag(tag);
}

void FEBioBoundarySection25::ParsePeriodicLinearConstraint(XMLTag& tag)
{
	FEModel* fem = GetFEModel();
	FEMesh& mesh = fem->GetMesh();
	FEPeriodicLinearConstraint plc(fem);

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

				FESurface* surf2 = fecore_alloc(FESurface, fem); feb->BuildSurface(*surf2, *spair->GetSecondarySurface());
				FESurface* surf1 = fecore_alloc(FESurface, fem); feb->BuildSurface(*surf1, *spair->GetPrimarySurface());
				plc.AddNodeSetPair(surf2->GetNodeList(), surf1->GetNodeList());
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
	assert(false);
	throw XMLReader::InvalidTag(tag);

/*	FEModel* fem = GetFEModel();
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

				FESurface* surf2 = fecore_alloc(FESurface, fem); feb->BuildSurface(*surf2, *spair->GetSecondarySurface());
				FESurface* surf1 = fecore_alloc(FESurface, fem); feb->BuildSurface(*surf1, *spair->GetPrimarySurface());
				plc.AddNodeSetPair(surf2->GetNodeList(), surf1->GetNodeList());
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
	*/
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

				FESurface& s = *(ntype == 1? pci->GetSecondarySurface() : pci->GetPrimarySurface());
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

	FEModelBuilder* feb = GetBuilder();

	// get the type attribute
	const char* szt = tag.AttributeValue("type");

	// Not all contact interfaces can be parsed automatically.
	// First, check all these special cases.
	if (strcmp(szt, "rigid_wall") == 0)
	{
		// --- R I G I D   W A L L   I N T E R F A C E ---

		FESurfacePairConstraint* ps = fecore_new<FESurfacePairConstraint>(szt, GetFEModel());
		if (ps)
		{
			fem.AddSurfacePairConstraint(ps);

			++tag;
			do
			{
				if (ReadParameter(tag, ps) == false)
				{
					if (tag == "surface")
					{
						FESurface& s = *ps->GetPrimarySurface();

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
		FENodalBC* prn = 0;
		FENodeSet* ns = 0;
		for (int i=0; i<nrn; ++i)
		{
			id = atoi(tag.AttributeValue("id"))-1;
			rb = atoi(tag.AttributeValue("rb"));

			if ((prn == 0) || (rb != rbp))
			{
				prn = fecore_new_class<FENodalBC>("FERigidNodeSet", &fem);

				prn->SetParameter("rb", rb);

				ns = new FENodeSet(&fem);
				prn->SetNodeSet(ns);

				// the default shell bc depends on the shell formulation
				// hinged shell = 0
				// clamped shell = 1
				prn->SetParameter("clamp_shells", feb->m_default_shell == OLD_SHELL ? 0 : 1);

				feb->AddBC(prn);
				rbp = rb;
			}
			ns->Add(id);

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
				FEAugLagLinearConstraint* pLC = fecore_alloc(FEAugLagLinearConstraint, &fem);

				++tag;
				do
				{
					int node, bc;
					double val;
					if (tag == "node")
					{
						tag.value(val);

						const char* szid = tag.AttributeValue("id");
						node = atoi(szid);

						const char* szbc = tag.AttributeValue("bc");
						if      (strcmp(szbc, "x") == 0) bc = 0;
						else if (strcmp(szbc, "y") == 0) bc = 1;
						else if (strcmp(szbc, "z") == 0) bc = 2;
						else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

						pLC->AddDOF(node, bc, val);
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
		FESurfacePairConstraint* pci = fecore_new<FESurfacePairConstraint>(szt, GetFEModel());
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
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	FEModelBuilder* feb = GetBuilder();
	int NMAT = fem.Materials();

	// get the rigid body material ID
	int rb = -1;
	tag.AttributeValue("rb", rb);

	// make sure we have a valid rigid body reference
	if ((rb <= 0)||(rb>NMAT)) throw XMLReader::InvalidAttributeValue(tag, "rb", tag.AttributeValue("rb"));

	// get the nodeset
	const char* szset = tag.AttributeValue("node_set");
	FENodeSet* nodeSet = mesh.FindNodeSet(szset);
	if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

	// create new rigid node set
	FENodalBC* prn = fecore_new_class<FENodalBC>("FERigidNodeSet", &fem);

	// the default shell bc depends on the shell formulation
	prn->SetParameter("clamped_shells", feb->m_default_shell == OLD_SHELL ? 0 : 1);
	prn->SetParameter("rb", rb);
	prn->SetNodeSet(nodeSet);

	// add it to the model
	feb->AddBC(prn);
	
	// read the parameter list
	ReadParameterList(tag, prn);
}

//-----------------------------------------------------------------------------
// The rigid body "constraints" are moved to the Boundary section in 2.5
void FEBioBoundarySection25::ParseRigidBody(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEModelBuilder& feb = *GetBuilder();

	const char* szm = tag.AttributeValue("mat");
	assert(szm);

	// get the material ID
	int nmat = atoi(szm);
	if ((nmat <= 0) || (nmat > fem.Materials())) throw XMLReader::InvalidAttributeValue(tag, "mat", szm);

	if (tag.isleaf()) return;

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

			double val = 0.0;
			value(tag, val);

			// create the rigid displacement constraint
			FEStepComponent* pDC = fecore_new_class<FEBoundaryCondition>("FERigidPrescribedOld", &fem);
			feb.AddRigidComponent(pDC);

			pDC->SetParameter("rb", nmat);
			pDC->SetParameter("dof", bc);
			pDC->SetParameter("relative", brel);
			pDC->SetParameter("value", val);

			// assign a load curve
			if (lc >= 0)
			{
				FEParam* p = pDC->GetParameter("value");
				if (p == nullptr) throw XMLReader::InvalidTag(tag);
				GetFEModel()->AttachLoadController(p, lc);
			}
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
			int ntype = 0; // FERigidBodyForce::FORCE_LOAD;
			const char* sztype = tag.AttributeValue("type", true);
			if (sztype)
			{
				if      (strcmp(sztype, "ramp"  ) == 0) ntype = 2; //FERigidBodyForce::FORCE_TARGET;
				else if (strcmp(sztype, "follow") == 0) ntype = 1; //FERigidBodyForce::FORCE_FOLLOW;
				else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
			}

			// get the loadcurve
			const char* szlc = tag.AttributeValue("lc", true);
			int lc = -1;
			if (szlc) lc = atoi(szlc) - 1;

			// make sure there is a loadcurve for type=0 forces
			if ((ntype == 0)&&(lc==-1)) throw XMLReader::MissingAttribute(tag, "lc");

			double val = 0.0;
			value(tag, val);

			// create the rigid body force/moment
			FEModelLoad* pFC = nullptr;
			if (bc < 3)
			{
				pFC = fecore_new<FEModelLoad>("rigid_force", &fem);

				pFC->SetParameter("load_type", ntype);
				pFC->SetParameter("rb", nmat);
				pFC->SetParameter("dof", bc);
				pFC->SetParameter("value", val);
			}
			else
			{
				pFC = fecore_new<FEModelLoad>("rigid_moment", &fem);

				pFC->SetParameter("rb", nmat);
				pFC->SetParameter("dof", bc - 3);
				pFC->SetParameter("value", val);
			}

			feb.AddModelLoad(pFC);

			if (lc >= 0)
			{
				FEParam* p = pFC->GetParameter("value");
				if (p == nullptr) throw XMLReader::InvalidTag(tag);
				GetFEModel()->AttachLoadController(p, lc);
			}
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
			FEBoundaryCondition* pBC = fecore_new_class<FEBoundaryCondition>("FERigidFixedBCOld",  &fem);
			feb.AddRigidComponent(pBC);

			pBC->SetParameter("rb", nmat);

			vector<int> dofs; dofs.push_back(bc);
			pBC->SetParameter("dofs", dofs);

		}
		else if (tag == "initial_velocity")
		{
			// get the initial velocity
			vec3d v;
			value(tag, v);

			// create the initial condition
			FEStepComponent* pic = fecore_new_class<FEInitialCondition>("FERigidBodyVelocity", &fem);
			pic->SetParameter("rb", nmat);
			pic->SetParameter("value", v);

			// add to model
			feb.AddRigidComponent(pic);
		}
		else if (tag == "initial_angular_velocity")
		{
			// get the initial angular velocity
			vec3d w;
			value(tag, w);

			// create the initial condition
			FEStepComponent* pic = fecore_new_class<FEInitialCondition>("FERigidBodyAngularVelocity", &fem);
			pic->SetParameter("rb", nmat);
			pic->SetParameter("value", w);

			// add to model
			feb.AddRigidComponent(pic);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
