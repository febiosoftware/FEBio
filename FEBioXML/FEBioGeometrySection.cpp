#include "stdafx.h"
#include "FEBioGeometrySection.h"
#include "FECore/FESolidDomain.h"
#include "FECore/FEShellDomain.h"
#include "FECore/FETrussDomain.h"
#include "FECore/FEDomain2D.h"
#include "FECore/FEModel.h"
#include "FEBioMech/FEElasticMaterial.h"
#include "FECore/FECoreKernel.h"
#include "FEBioMech/FEElasticMixture.h"

//-----------------------------------------------------------------------------
//!  Parses the geometry section from the xml file
//!
void FEBioGeometrySection::Parse(XMLTag& tag)
{
	m_pim->m_maxid = 0;
	int nversion = m_pim->Version();

	if (nversion < 0x0200)
	{
		++tag;
		do
		{
			if      (tag == "Nodes"      ) ParseNodeSection       (tag);
			else if (tag == "Elements"   ) ParseElementSection    (tag);
			else if (tag == "ElementData") ParseElementDataSection(tag);
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}
	else
	{
		++tag;
		do
		{
			if      (tag == "Nodes"      ) ParseNodeSection       (tag);
			else if (tag == "Elements"   ) ParseElementSection20  (tag);
			else if (tag == "NodeSet"    ) ParseNodeSetSection    (tag);
			else if (tag == "Surface"    ) ParseSurfaceSection    (tag);
			else if (tag == "Edge"       ) ParseEdgeSection       (tag);
			else if (tag == "ElementSet" ) ParseElementSetSection (tag);
			else if (tag == "DiscreteSet") 
			{
				if (nversion >= 0x0205) ParseDiscreteSetSection(tag);
				else throw XMLReader::InvalidTag(tag);
			}
			else if (tag == "SurfacePair")
			{
				if (nversion >= 0x0205) ParseSurfacePairSection(tag);
				else throw XMLReader::InvalidTag(tag);
			}
			else if (tag == "NodeSetPair")
			{
				if (nversion >= 0x0205) ParseNodeSetPairSection(tag);
				else throw XMLReader::InvalidTag(tag);
			}
			else if (tag == "ElementData")
			{
				// This section is no longer supported since 2.5 and replaced
				// with the top level MeshData section.
				if (nversion < 0x0205) ParseElementDataSection(tag);
				else throw XMLReader::InvalidTag(tag);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}

	// At this point the mesh is completely read in.
	// Now we can allocate the degrees of freedom.
	// NOTE: We do this here since the mesh no longer automatically allocates the dofs.
	//       At some point I want to be able to read the mesh before deciding any physics.
	//       When that happens I'll have to move this elsewhere.
	FEModel& fem = *GetFEModel();
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	fem.GetMesh().SetDOFS(MAX_DOFS);
}

//-----------------------------------------------------------------------------
//! Reads the Nodes section of the FEBio input file
void FEBioGeometrySection::ParseNodeSection(XMLTag& tag)
{
	FEMesh& mesh = *m_pim->GetFEMesh();
	int N0 = mesh.Nodes();

	// get the largest nodal ID
	// (It is assumed that nodes are sorted by ID so the last node should have the largest ID)
	int max_id = 0;
	if (N0 > 0) max_id = mesh.Node(N0-1).GetID();

	// first we need to figure out how many nodes there are
	XMLTag t(tag);
	int nodes = tag.children();

	// see if this list defines a set
	const char* szref = (m_pim->Version() >= 0x0205 ? "name" : "set");
	const char* szl = tag.AttributeValue(szref, true);
	FENodeSet* ps = 0;
	if (szl)
	{
		ps = new FENodeSet(&mesh);
		ps->SetName(szl);
		ps->create(nodes);
		mesh.AddNodeSet(ps);
	}

	// resize node's array
	mesh.AddNodes(nodes);

	// read nodal coordinates
	++tag;
	for (int i=0; i<nodes; ++i)
	{
		FENode& node = mesh.Node(N0 + i);
		m_pim->value(tag, node.m_r0);
		node.m_rt = node.m_r0;

		// get the nodal ID
		int nid = -1;
		tag.AttributeValue("id", nid);

		// Make sure it is valid
		if (nid <= max_id) throw XMLReader::InvalidAttributeValue(tag, "id");

		// set the ID
		node.SetID(nid);
		max_id = nid;

		// go on to the next node
		++tag;
	}

	// If a node-set is defined add these nodes to the node-set
	if (ps)
	{
		for (int i=0; i<nodes; ++i) (*ps)[i] = N0+i;
	}

	// tell the file reader to rebuild the node ID table
	m_pim->BuildNodeList();
}

//-----------------------------------------------------------------------------
//! Get the element type from a XML tag
FE_Element_Spec FEBioGeometrySection::ElementSpec(const char* sztype)
{
	// determine the element shape 
	FE_Element_Shape eshape = FE_ELEM_INVALID_SHAPE;
	if      (strcmp(sztype, "hex8"  ) == 0) eshape = ET_HEX8;
	else if (strcmp(sztype, "hex20" ) == 0) eshape = ET_HEX20;
	else if (strcmp(sztype, "hex27" ) == 0) eshape = ET_HEX27;
	else if (strcmp(sztype, "penta6") == 0) eshape = ET_PENTA6;
	else if (strcmp(sztype, "tet4"  ) == 0) eshape = ET_TET4;
	else if (strcmp(sztype, "tet10" ) == 0) eshape = ET_TET10;
	else if (strcmp(sztype, "tet15" ) == 0) eshape = ET_TET15;
	else if (strcmp(sztype, "quad4" ) == 0) eshape = ET_QUAD4;
	else if (strcmp(sztype, "quad8" ) == 0) eshape = ET_QUAD8;
	else if (strcmp(sztype, "quad9" ) == 0) eshape = ET_QUAD9;
	else if (strcmp(sztype, "tri3"  ) == 0) eshape = ET_TRI3;
	else if (strcmp(sztype, "tri6"  ) == 0) eshape = ET_TRI6;
	else if (strcmp(sztype, "truss2") == 0) eshape = ET_TRUSS2;
	else
	{
		// new way for defining element type and integration rule at the same time
		// this is useful for multi-step analyses where the geometry is read in before the control section.
		if      (strcmp(sztype, "TET10G4"     ) == 0) { eshape = ET_TET10; m_pim->m_ntet10 = FE_TET10G4; }
		else if (strcmp(sztype, "TET10G8"     ) == 0) { eshape = ET_TET10; m_pim->m_ntet10 = FE_TET10G8; }
		else if (strcmp(sztype, "TET10GL11"   ) == 0) { eshape = ET_TET10; m_pim->m_ntet10 = FE_TET10GL11; }
		else if (strcmp(sztype, "TET10G4_S3"  ) == 0) { eshape = ET_TET10; m_pim->m_ntet10 = FE_TET10G4; m_pim->m_ntri6 = FE_TRI6G3; }
		else if (strcmp(sztype, "TET10G8_S3"  ) == 0) { eshape = ET_TET10; m_pim->m_ntet10 = FE_TET10G8; m_pim->m_ntri6 = FE_TRI6G3; }
		else if (strcmp(sztype, "TET10GL11_S3") == 0) { eshape = ET_TET10; m_pim->m_ntet10 = FE_TET10GL11; m_pim->m_ntri6 = FE_TRI6G3; }
		else if (strcmp(sztype, "TET10G4_S4"  ) == 0) { eshape = ET_TET10; m_pim->m_ntet10 = FE_TET10G4; m_pim->m_ntri6 = FE_TRI6G4; }
		else if (strcmp(sztype, "TET10G8_S4"  ) == 0) { eshape = ET_TET10; m_pim->m_ntet10 = FE_TET10G8; m_pim->m_ntri6 = FE_TRI6G4; }
		else if (strcmp(sztype, "TET10GL11_S4") == 0) { eshape = ET_TET10; m_pim->m_ntet10 = FE_TET10GL11; m_pim->m_ntri6 = FE_TRI6G4; }
		else if (strcmp(sztype, "TET10G4_S7"  ) == 0) { eshape = ET_TET10; m_pim->m_ntet10 = FE_TET10G4; m_pim->m_ntri6 = FE_TRI6G7; }
		else if (strcmp(sztype, "TET10G8_S7"  ) == 0) { eshape = ET_TET10; m_pim->m_ntet10 = FE_TET10G8; m_pim->m_ntri6 = FE_TRI6G7; }
		else if (strcmp(sztype, "TET10GL11_S7") == 0) { eshape = ET_TET10; m_pim->m_ntet10 = FE_TET10GL11; m_pim->m_ntri6 = FE_TRI6G7; }
		else if (strcmp(sztype, "TET15G8"     ) == 0) { eshape = ET_TET15; m_pim->m_ntet15 = FE_TET15G8; }
		else if (strcmp(sztype, "TET15G11"    ) == 0) { eshape = ET_TET15; m_pim->m_ntet15 = FE_TET15G11; }
		else if (strcmp(sztype, "TET15G15"    ) == 0) { eshape = ET_TET15; m_pim->m_ntet15 = FE_TET15G15; }
		else if (strcmp(sztype, "TET15G8_S3"  ) == 0) { eshape = ET_TET15; m_pim->m_ntet15 = FE_TET15G8; m_pim->m_ntri7 = FE_TRI7G3; }
		else if (strcmp(sztype, "TET15G11_S3" ) == 0) { eshape = ET_TET15; m_pim->m_ntet15 = FE_TET15G11; m_pim->m_ntri7 = FE_TRI7G3; }
		else if (strcmp(sztype, "TET15G15_S3" ) == 0) { eshape = ET_TET15; m_pim->m_ntet15 = FE_TET15G15; m_pim->m_ntri7 = FE_TRI7G3; }
		else if (strcmp(sztype, "TET15G8_S4"  ) == 0) { eshape = ET_TET15; m_pim->m_ntet15 = FE_TET15G8; m_pim->m_ntri7 = FE_TRI7G4; }
		else if (strcmp(sztype, "TET15G11_S4" ) == 0) { eshape = ET_TET15; m_pim->m_ntet15 = FE_TET15G11; m_pim->m_ntri7 = FE_TRI7G4; }
		else if (strcmp(sztype, "TET15G15_S4" ) == 0) { eshape = ET_TET15; m_pim->m_ntet15 = FE_TET15G15; m_pim->m_ntri7 = FE_TRI7G4; }
		else if (strcmp(sztype, "TET15G8_S7"  ) == 0) { eshape = ET_TET15; m_pim->m_ntet15 = FE_TET15G8; m_pim->m_ntri7 = FE_TRI7G7; }
		else if (strcmp(sztype, "TET15G11_S7" ) == 0) { eshape = ET_TET15; m_pim->m_ntet15 = FE_TET15G11; m_pim->m_ntri7 = FE_TRI7G7; }
		else if (strcmp(sztype, "TET15G15_S7" ) == 0) { eshape = ET_TET15; m_pim->m_ntet15 = FE_TET15G15; m_pim->m_ntri7 = FE_TRI7G7; }
		else
		{
			throw FEBioImport::InvalidElementType();
		}
	}

	// this is a hack to choose between 2D elements and shell elements.
	// NOTE: This is only used by quad/tri elements.
	// TODO: find a better way
	int NDIM = 3;
	if (strcmp(m_pim->m_szmod, "fluid") == 0) NDIM = 2;

	// determine the element type
	FE_Element_Type etype = FE_ELEM_INVALID_TYPE;
	switch (eshape)
	{
	case ET_HEX8  : etype = m_pim->m_nhex8; break;
	case ET_PENTA6: etype = FE_PENTA6G6; break;
	case ET_TET4  : etype = m_pim->m_ntet4; break;
	case ET_TET10 : etype = m_pim->m_ntet10; break;
	case ET_TET15 : etype = m_pim->m_ntet15; break;
	case ET_HEX20 : etype = FE_HEX20G27; break;
	case ET_HEX27 : etype = FE_HEX27G27; break;
	case ET_QUAD4 : etype = (NDIM == 3 ? FE_SHELL_QUAD : FE2D_QUAD4G4); break;
	case ET_TRI3  : etype = (NDIM == 3 ? FE_SHELL_TRI  : FE2D_TRI3G1 ); break;
    case ET_TRI6  : etype = (NDIM == 3 ? FE_SHELL_TRI6 : FE2D_TRI6G3); break;
	case ET_QUAD8 : etype = (NDIM == 3 ? FE_SHELL_QUAD8 : FE2D_QUAD8G9); break;
	case ET_QUAD9 : etype = FE2D_QUAD9G9; break;
	case ET_TRUSS2: etype = FE_TRUSS; break;
	default:
		throw FEBioImport::InvalidElementType();
	}

	// determine the element class
	FE_Element_Class eclass = FEElementLibrary::GetElementClass(etype);

	// return the spec
	FE_Element_Spec spec;
	spec.eclass = eclass;
	spec.eshape = eshape;
	spec.etype  = etype;
	spec.m_bthree_field_hex = m_pim->m_b3field_hex;
	spec.m_bthree_field_tet = m_pim->m_b3field_tet;
	spec.m_but4 = m_pim->m_but4;

	// Make sure this is a valid element specification
	if (FEElementLibrary::IsValid(spec) == false) FEBioImport::InvalidElementType();

	return spec;
}

//-----------------------------------------------------------------------------
//! find the domain type for the element and material type
FEDomain* FEBioGeometrySection::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	// get the domain type
	FECoreKernel& febio = FECoreKernel::GetInstance();
	return febio.CreateDomain(spec, pm, pmat);
}

//-----------------------------------------------------------------------------
//! This function reads the Element section from the FEBio input file. It also
//! creates the domain classes which store the element data. A domain is defined
//! by the module (structural, poro, heat, etc), the element type (solid, shell,
//! etc.) and the material. 
//!
void FEBioGeometrySection::ParseElementSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = *m_pim->GetFEMesh();

	// first we need to figure out how many elements 
	// and how many domains there are
	vector<FEDOMAIN> dom;
	vector<int>	ED; ED.reserve(1000);
	XMLTag t(tag); ++t;
	int elems = 0, i;
	while (!t.isend())
	{
		// get the material ID
		const char* szmat = t.AttributeValue("mat");
		int nmat = atoi(szmat)-1;
		if ((nmat < 0) || (nmat >= fem.Materials())) throw FEBioImport::InvalidMaterial(elems+1);

		// get the element type
		FE_Element_Spec spec = ElementSpec(t.Name());

		// get the element ID
		int nid = -1;
		t.AttributeValue("id", nid);

		// keep track of the largest element ID
		if (nid > m_pim->m_maxid) m_pim->m_maxid = nid;

		// find a domain for this element
		int ndom = -1;
		for (i=0; i<(int) dom.size(); ++i)
		{
			FEDOMAIN& d = dom[i];
			if ((d.mat == nmat) && (d.elem.eshape == spec.eshape))
			{
				ndom = i;
				d.nel++;
				break;
			}
		}
		if (ndom == -1)
		{
			FEDOMAIN d;
			d.mat = nmat;
			d.elem = spec;
			d.nel = 1;
			ndom = (int) dom.size();
			dom.push_back(d);
		}

		ED.push_back(ndom);
		elems++;
		++t;
	}

	// create the domains
	for (i=0; i<(int) dom.size(); ++i)
	{
		FEDOMAIN& d = dom[i];

		// get material class
		FEMaterial* pmat = fem.GetMaterial(d.mat);

		// create the new domain
		FEDomain* pdom = CreateDomain(d.elem, &mesh, pmat);
		if (pdom == 0) throw FEBioImport::FailedCreatingDomain();

		// add it to the mesh
		assert(d.nel);
		pdom->Create(d.nel, d.elem.etype);
		mesh.AddDomain(pdom);

		// we reset the nr of elements since we'll be using 
		// that variable is a counter in the next loop
		d.nel = 0;
	}

	// read element data
	++tag;
	int nid = 1;
	for (i=0; i<elems; ++i, ++nid)
	{
		int nd = ED[i];
		int ne = dom[nd].nel++;

		// get the domain to which this element belongs
		FEDomain& domi = mesh.Domain(nd);

		// get the material ID
		int nmat = atoi(tag.AttributeValue("mat"))-1;
		domi.SetMatID(nmat);

		// get the material class
		FEMaterial* pmat = fem.GetMaterial(nmat);
		assert(pmat == domi.GetMaterial());

		// determine element shape
		FE_Element_Spec espec = ElementSpec(tag.Name());
		if (espec.etype != dom[ED[i]].elem.etype) throw XMLReader::InvalidTag(tag);
		ReadElement(tag, domi.ElementRef(ne), nid);

		// go to next tag
		++tag;
	}

	// assign material point data
	for (i=0; i<mesh.Domains(); ++i)
	{
		FEDomain& d = mesh.Domain(i);
		d.CreateMaterialPointData();
	}
}

//-----------------------------------------------------------------------------
//! This function reads the Element section from the FEBio input file. It also
//! creates the domain classes which store the element data. A domain is defined
//! by the module (structural, poro, heat, etc), the element type (solid, shell,
//! etc.) and the material. 
//!
void FEBioGeometrySection::ParseElementSection20(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	int NDOM = (int) m_dom.size();

	// get the material ID
	const char* szmat = tag.AttributeValue("mat");
	int nmat = atoi(szmat)-1;
	if ((nmat < 0) || (nmat >= fem.Materials())) throw FEBioImport::InvalidDomainMaterial();

	// get the name
	const char* szref = (m_pim->Version() >= 0x0205 ? "name" : "elset");
	const char* szname = tag.AttributeValue(szref, true);

	// get the element type
	const char* sztype = tag.AttributeValue("type");
	FE_Element_Spec espec = ElementSpec(sztype);

	// get the domain's material class
	FEMaterial* pmat = fem.GetMaterial(nmat);

	// create the new domain
	FEDomain* pdom = CreateDomain(espec, &mesh, pmat);
	if (pdom == 0) throw FEBioImport::FailedCreatingDomain();
	FEDomain& dom = *pdom;
	dom.SetName(szname);

	// count elements
	int elems = tag.children();
	assert(elems);

	// add domain it to the mesh
	pdom->Create(elems, espec.etype);
	pdom->SetMatID(nmat);
	mesh.AddDomain(pdom);
	int nd = NDOM;

	// for named domains, we'll also create an element set
	FEElementSet* pg = 0;
	if (szname)
	{
		pg = new FEElementSet(&mesh);
		pg->SetName(szname);
		pg->create(elems);
		mesh.AddElementSet(pg);
	}

	// read element data
	++tag;
	for (int i=0; i<elems; ++i)
	{
		if ((tag == "elem")==false) throw XMLReader::InvalidTag(tag);

		// get the element ID
		int nid;
		tag.AttributeValue("id", nid);

		// Make sure element IDs increase
//		if (nid <= m_pim->m_maxid) throw XMLReader::InvalidAttributeValue(tag, "id");

		// keep track of the largest element ID
		// (which by assumption is the ID that was just read in)
		m_pim->m_maxid = nid;

		// add to the element set (if we have one)
		if (pg) (*pg)[i] = nid;

		// read the element data
		ReadElement(tag, dom.ElementRef(i), nid);

		// go to next tag
		++tag;
	}

	// assign material point data
	dom.CreateMaterialPointData();
}

//-----------------------------------------------------------------------------
//! This function reads the Element section from the FEBio input file. It also
//! creates the domain classes which store the element data. A domain is defined
//! by the module (structural, poro, heat, etc), the element type (solid, shell,
//! etc.) and the material. 
//!
void FEBioGeometrySection::ParseMesh(XMLTag& tag)
{
	FEMesh& mesh = *m_pim->GetFEMesh();

	// first we need to figure out how many elements 
	// and how many domains there are
	vector<FEDOMAIN> dom;
	vector<int>	ED; ED.reserve(1000);
	XMLTag t(tag); ++t;
	int elems = 0, i;
	while (!t.isend())
	{
		// get the element type
		FE_Element_Spec espec = ElementSpec(t.Name());

		// find a domain for this element
		int ndom = -1;
		for (i=0; i<(int) dom.size(); ++i)
		{
			FEDOMAIN& d = dom[i];
			if (d.elem == espec)
			{
				ndom = i;
				d.nel++;
				break;
			}
		}
		if (ndom == -1)
		{
			FEDOMAIN d;
			d.mat = 0;
			d.elem = espec;
			d.nel = 1;
			ndom = (int) dom.size();
			dom.push_back(d);
		}

		ED.push_back(ndom);
		elems++;
		++t;
	}

	// create the domains
	for (i=0; i<(int) dom.size(); ++i)
	{
		FEDOMAIN& d = dom[i];

		// create the new domain
		FEDomain* pdom = CreateDomain(d.elem, &mesh, 0);
		if (pdom == 0) throw FEBioImport::FailedCreatingDomain();


		// add it to the mesh
		assert(d.nel);
		pdom->Create(d.nel, d.elem.etype);
		pdom->SetMatID(0);
		mesh.AddDomain(pdom);

		// we reset the nr of elements since we'll be using 
		// that variable is a counter in the next loop
		d.nel = 0;
	}

	// read element data
	++tag;
	int nid = 1;
	for (i=0; i<elems; ++i, ++nid)
	{
		int nd = ED[i];
		int ne = dom[nd].nel++;

		// get the domain to which this element belongs
		FEDomain& dom = mesh.Domain(nd);

		// determine element type
		FE_Element_Spec espec = ElementSpec(tag.Name());
		ReadElement(tag, dom.ElementRef(ne), nid); break;

		// go to next tag
		++tag;
	}
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection::ReadElement(XMLTag &tag, FEElement& el, int nid)
{
	el.SetID(nid);
	int n[FEElement::MAX_NODES];
	tag.value(n,el.Nodes());
	m_pim->GlobalToLocalID(n, el.Nodes(), el.m_node);
}

//-----------------------------------------------------------------------------
void set_element_fiber(FEElement& el, const vec3d& v)
{
	// normalize fiber
	vec3d a = v;
	a.unit();

	// set up a orthonormal coordinate system
	vec3d b(0,1,0);
	if (fabs(fabs(a*b) - 1) < 1e-7) b = vec3d(0,0,1);
	vec3d c = a^b;
	b = c^a;

	// make sure they are unit vectors
	b.unit();
	c.unit();

	for (int i=0; i<el.GaussPoints(); ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3d& m = pt.m_Q;
		m.zero();
		m[0][0] = a.x; m[0][1] = b.x; m[0][2] = c.x;
		m[1][0] = a.y; m[1][1] = b.y; m[1][2] = c.y;
		m[2][0] = a.z; m[2][1] = b.z; m[2][2] = c.z;
	}
}

//-----------------------------------------------------------------------------
void set_element_mat_axis(FEElement& el, const vec3d& v1, const vec3d& v2)
{
	vec3d a = v1;
	vec3d d = v2;
				
	vec3d c = a^d;
	vec3d b = c^a;

	// normalize
	a.unit();
	b.unit();
	c.unit();

	// assign to element
	for (int i=0; i<el.GaussPoints(); ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		pt.m_Q = mat3d(a, b, c);
	}
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection::ParseElementData(FEElement& el, XMLTag& tag)
{
	vec3d a, d;
	if (tag == "fiber")
	{
		// read the fiber direction
		m_pim->value(tag, a);
		set_element_fiber(el, a);
	}
	else if (tag == "mat_axis")
	{
		++tag;
		do
		{
			if      (tag == "a") m_pim->value(tag, a);
			else if (tag == "d") m_pim->value(tag, d);
			else throw XMLReader::InvalidTag(tag);

			++tag;
		}
		while (!tag.isend());
		set_element_mat_axis(el, a, d);
	}
	else if (tag == "thickness")
	{
		if (el.Class() != FE_ELEM_SHELL) throw XMLReader::InvalidTag(tag);
		FEShellElement& shell = static_cast<FEShellElement&> (el);

		// read shell thickness
		tag.value(&shell.m_h0[0],shell.Nodes());
	}
	else if (tag == "area")
	{
		if (el.Class() != FE_ELEM_TRUSS) throw XMLReader::InvalidTag(tag);
		FETrussElement& truss = static_cast<FETrussElement&>(el);

		// read truss area
		m_pim->value(tag, truss.m_a0);
	}
	else
	{
		for (int i=0; i<el.GaussPoints(); ++i)
		{
			FEMaterialPoint* pt = el.GetMaterialPoint(i);
			while (pt)
			{
				FEParameterList& pl = pt->GetParameterList();
				if (m_pim->ReadParameter(tag, pl)) break;

				FEElasticMixtureMaterialPoint* mPt = dynamic_cast<FEElasticMixtureMaterialPoint*>(pt);

				bool tagFound=false;
				if(mPt)
				{
					vector<FEMaterialPoint*> mPtV = mPt->m_mp;
					for (int i=0; i<(int)mPtV.size(); ++i)
					{
						FEParameterList& pl = mPtV[i]->GetParameterList();
						if (m_pim->ReadParameter(tag, pl))
						{
							tagFound=true;
							break;
						}
					}
				}

				if(tagFound) break;

				pt = pt->Next();
				if (pt == 0) throw XMLReader::InvalidTag(tag);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Reads the ElementData section from the FEBio input file

void FEBioGeometrySection::ParseElementDataSection(XMLTag& tag)
{
	int i;

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = *m_pim->GetFEMesh();

	// get the total nr of elements
	int nelems = mesh.Elements();

	//make sure we've read the element section
	if (nelems == 0) throw XMLReader::InvalidTag(tag);

	// create the pelem array
	vector<FEElement*> pelem;
	pelem.assign(nelems, static_cast<FEElement*>(0));

	for (int nd=0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& d = mesh.Domain(nd);
		for (i=0; i<d.Elements(); ++i)
		{
			FEElement& el = d.ElementRef(i);
			assert(pelem[el.GetID()-1] == 0);
			pelem[el.GetID()-1] = &el;
		}
	}

	// read additional element data
	++tag;
	do
	{
		// make sure this is an "element" tag
		if (tag == "element")
		{
			// get the element number
			const char* szid = tag.AttributeValue("id");
			int n = atoi(szid)-1;

			// make sure the number is valid
			if ((n<0) || (n>=nelems)) throw XMLReader::InvalidAttributeValue(tag, "id", szid);

			// get a pointer to the element
			FEElement* pe = pelem[n];

			++tag;
			do
			{
				ParseElementData(*pe, tag);
				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "elset")
		{
			const char* szname = tag.AttributeValue("set");
			// find domain with this name
			FEElementSet* pset = mesh.FindElementSet(szname);
			if (pset == 0) throw XMLReader::InvalidAttributeValue(tag, "set", szname);

			++tag;
			do
			{
				int n = pset->size();
				for (int i=0; i<n; ++i)
				{
					// get a pointer to the element
					int nid = (*pset)[i] - 1;
					FEElement* pe = pelem[nid];
					ParseElementData(*pe, tag);
				}
				++tag;
			}
			while (!tag.isend());
		}
		else throw XMLReader::InvalidTag(tag);
		
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
//! Reads the Geometry::Groups section of the FEBio input file

void FEBioGeometrySection::ParseNodeSetSection(XMLTag& tag)
{
	int nversion = m_pim->Version();
	const char* szatt = (nversion < 0x0205 ? "set" : "node_set");

	// read the node set
	FENodeSet* pns = m_pim->ParseNodeSet(tag, szatt);
	if (pns == 0) throw XMLReader::InvalidTag(tag);
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection::ParseDiscreteSetSection(XMLTag& tag)
{
	// get the mesh
	FEMesh& mesh = *m_pim->GetFEMesh();

	// get the name
	const char* szname = tag.AttributeValue("name");

	// create the discrete element set
	FEDiscreteSet* ps = new FEDiscreteSet(&mesh);
	ps->SetName(szname);
	mesh.AddDiscreteSet(ps);

	// read the node pairs
	++tag;
	do
	{
		if (tag == "delem")
		{
			int n[2];
			tag.value(n, 2);
			n[0] -= 1; n[1] -= 1;
			ps->add(n[0], n[1]);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
//! Reads a Geometry\Edge section.
void FEBioGeometrySection::ParseEdgeSection(XMLTag& tag)
{
	// get the mesh
	FEMesh& mesh = *m_pim->GetFEMesh();

	// get the number of nodes
	// (we use this for checking the node indices of the facets)
	int NN = mesh.Nodes();

	// get the required name attribute
	const char* szname = tag.AttributeValue("name");

	// count nr of segments
	int nsegs = tag.children();

	// allocate storage for segments
	FESegmentSet* ps = new FESegmentSet;
	ps->Create(nsegs);
	ps->SetName(szname);

	// add it to the mesh
	mesh.AddSegmentSet(ps);

	// read segments
	++tag;
	int nf[FEElement::MAX_NODES];
	for (int i=0; i<nsegs; ++i)
	{
		FESegmentSet::SEGMENT& line = ps->Segment(i);

		// set the facet type
		if      (tag == "line2") line.ntype = 2;
		else throw XMLReader::InvalidTag(tag);

		// we assume that the segment type also defines the number of nodes
		int N = line.ntype;
		tag.value(nf, N);
		for (int j=0; j<N; ++j) 
		{
			int nid = nf[j]-1;
			if ((nid<0)||(nid>= NN)) throw XMLReader::InvalidValue(tag);
			line.node[j] = nid;
		}

		++tag;
	}
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection::ParseSurfacePairSection(XMLTag& tag)
{
	FEBioImport::SurfacePair p;
	const char* szname = tag.AttributeValue("name");
	strcpy(p.szname, szname);

	FEMesh& mesh = *m_pim->GetFEMesh();

	++tag;
	do
	{
		if (tag == "master")
		{
			const char* sz = tag.AttributeValue("surface");
			p.pmaster = mesh.FindFacetSet(sz);
			if (p.pmaster == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);
		}
		else if (tag == "slave")
		{
			const char* sz = tag.AttributeValue("surface");
			p.pslave = mesh.FindFacetSet(sz);
			if (p.pslave == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	m_pim->AddSurfacePair(p);
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection::ParseNodeSetPairSection(XMLTag& tag)
{
	FEBioImport::NodeSetPair p;
	const char* szname = tag.AttributeValue("name");
	strcpy(p.szname, szname);

	FEMesh& mesh = *m_pim->GetFEMesh();

	++tag;
	do
	{
		if (tag == "master")
		{
			const char* sz = tag.AttributeValue("node_set");
			p.pmaster = mesh.FindNodeSet(sz);
			if (p.pmaster == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", sz);
		}
		else if (tag == "slave")
		{
			const char* sz = tag.AttributeValue("node_set");
			p.pslave = mesh.FindNodeSet(sz);
			if (p.pslave == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", sz);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());

	m_pim->AddNodeSetPair(p);
}

//-----------------------------------------------------------------------------
//! Reads a Geometry\Surface section.
void FEBioGeometrySection::ParseSurfaceSection(XMLTag& tag)
{
	// get the mesh
	FEMesh& mesh = *m_pim->GetFEMesh();

	// get the number of nodes
	// (we use this for checking the node indices of the facets)
	int NN = mesh.Nodes();

	// get the required name attribute
	const char* szname = tag.AttributeValue("name");

	// count nr of faces
	int faces = tag.children();

	// allocate storage for faces
	FEFacetSet* ps = new FEFacetSet;
	ps->Create(faces);
	ps->SetName(szname);

	// add it to the mesh
	mesh.AddFacetSet(ps);

	// read faces
	++tag;
	int nf[FEElement::MAX_NODES];
	for (int i=0; i<faces; ++i)
	{
		FEFacetSet::FACET& face = ps->Face(i);

		// set the facet type
		if      (tag == "quad4") face.ntype = 4;
		else if (tag == "tri3" ) face.ntype = 3;
		else if (tag == "tri6" ) face.ntype = 6;
		else if (tag == "tri7" ) face.ntype = 7;
		else if (tag == "quad8") face.ntype = 8;
		else if (tag == "quad9") face.ntype = 9;
		else throw XMLReader::InvalidTag(tag);

		// we assume that the facet type also defines the number of nodes
		int N = face.ntype;
		tag.value(nf, N);
		for (int j=0; j<N; ++j) 
		{
			int nid = nf[j]-1;
			if ((nid<0)||(nid>= NN)) throw XMLReader::InvalidValue(tag);
			face.node[j] = nid;
		}

		++tag;
	}
}

//-----------------------------------------------------------------------------
//! Reads the Geometry::Groups section of the FEBio input file

void FEBioGeometrySection::ParseElementSetSection(XMLTag& tag)
{
	FEMesh& mesh = *m_pim->GetFEMesh();

	// get the name attribute
	const char* szname = tag.AttributeValue("name");

	// create a new element set
	FEElementSet* pg = new FEElementSet(&mesh);
	pg->SetName(szname);

	vector<int> l;
	++tag;
	do
	{
		if (tag == "elem")
		{
			int nid = -1;
			tag.AttributeValue("id", nid);
			l.push_back(nid);
		}
		else { delete pg; throw XMLReader::InvalidTag(tag); }
		++tag;
	}
	while (!tag.isend());

	// only add non-empty element sets
	if (l.empty() == false)
	{
		// assign indices to element set
		int N = l.size();
		pg->create(N);
		for (int i=0; i<N; ++i) (*pg)[i] = l[i];

		// add the element set to the mesh
		mesh.AddElementSet(pg);
	}
	else delete pg;
}
