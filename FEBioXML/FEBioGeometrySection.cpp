#include "stdafx.h"
#include "FEBioGeometrySection.h"
#include "FECore/FESolidDomain.h"
#include "FECore/FEShellDomain.h"
#include "FECore/FETrussDomain.h"
#include "FECore/FEModel.h"
#include "FEBioMech/FEElasticMaterial.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
//!  Parses the geometry section from the xml file
//!
void FEBioGeometrySection::Parse(XMLTag& tag)
{
	m_pim->m_maxid = 0;

	if (m_pim->Version() < 0x0200)
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
			else if (tag == "ElementData") ParseElementDataSection(tag);
			else if (tag == "NodeSet"    ) ParseNodeSetSection    (tag);
			else if (tag == "Surface"    ) ParseSurfaceSection    (tag);
			else if (tag == "ElementSet" ) ParseElementSetSection (tag);
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}
}

//-----------------------------------------------------------------------------
//! Reads the Nodes section of the FEBio input file
void FEBioGeometrySection::ParseNodeSection(XMLTag& tag)
{
	FEMesh& mesh = *m_pim->GetFEMesh();
	int N0 = mesh.Nodes();

	// first we need to figure out how many nodes there are
	XMLTag t(tag);
	int nodes = 0;
	++t;
	while (!t.isend()) { nodes++; ++t; }

	// see if this list defines a set
	const char* szl = tag.AttributeValue("set", true);
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
		++tag;
	}

	// If a node-set is defined add these nodes to the node-set
	if (ps)
	{
		for (int i=0; i<nodes; ++i) (*ps)[i] = N0+i;
	}
}

//-----------------------------------------------------------------------------
//! Get the element type from a XML tag
FE_Element_Shape FEBioGeometrySection::ElementShape(XMLTag& t)
{
	if      (t=="hex8"  ) return ET_HEX8;
	else if (t=="hex20" ) return ET_HEX20;
	else if (t=="hex27" ) return ET_HEX27;
	else if (t=="penta6") return ET_PENTA6;
	else if (t=="tet4"  ) return ET_TET4;
	else if (t=="tet10" ) return ET_TET10;
	else if (t=="tet15" ) return ET_TET15;
	else if (t=="quad4" ) return ET_QUAD4;
	else if (t=="tri3"  ) return ET_TRI3;
	else if (t=="truss2") return ET_TRUSS2;
	else 
	{
		assert(false);
		throw XMLReader::InvalidTag(t);
	}
}

//-----------------------------------------------------------------------------
//! find the domain type for the element and material type
int FEBioGeometrySection::DomainType(FE_Element_Shape eshape, FEMaterial* pmat)
{
	// setup the element specs
	FE_Element_Spec spec;
	spec.eshape = eshape;
	spec.m_bthree_field = m_pim->m_b3field;
	spec.m_but4 = m_pim->m_but4;
	if      (eshape == ET_HEX8 ) spec.etype = m_pim->m_nhex8;
	else if (eshape == ET_TET4 ) spec.etype = m_pim->m_ntet4;
	else if (eshape == ET_TET10) spec.etype = m_pim->m_ntet10;
	else if (eshape == ET_TET15) spec.etype = m_pim->m_ntet15;
	else if (eshape == ET_TRI3 ) spec.etype = m_pim->m_ntri3;
	
	// get the domain type
	FECoreKernel& febio = FECoreKernel::GetInstance();
	return febio.GetDomainType(spec, pmat);
}

//-----------------------------------------------------------------------------
//! Create a particular type of domain
FEDomain* FEBioGeometrySection::CreateDomain(int ntype, FEMesh* pm, FEMaterial* pmat)
{
	// create a new domain based on the type
	FECoreKernel& febio = FECoreKernel::GetInstance();
	FEDomain* pd = febio.CreateDomain(ntype, pm, pmat);

	// return the domain
	return pd;
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
		FE_Element_Shape etype = ElementShape(t);

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
			if ((d.mat == nmat) && (d.elem == etype))
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
			d.elem = etype;
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

		// then, find the domain type depending on the 
		// element and material types
		int ntype = DomainType(d.elem, pmat);
		if (ntype == 0) throw FEBioImport::InvalidDomainType();

		// create the new domain
		FEDomain* pdom = CreateDomain(ntype, &mesh, pmat);
		if (pdom == 0) throw FEBioImport::FailedCreatingDomain();


		// add it to the mesh
		assert(d.nel);
		pdom->create(d.nel);
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

		// get the material ID
		int nmat = atoi(tag.AttributeValue("mat"))-1;

		// get the material class
		FEMaterial* pmat = fem.GetMaterial(nmat);
		assert(pmat == dom.GetMaterial());

		// determine element type
		int etype = -1;
		if      (tag == "hex8"  ) etype = ET_HEX8;
		else if (tag == "hex20" ) etype = ET_HEX20;
		else if (tag == "hex27" ) etype = ET_HEX27;
		else if (tag == "penta6") etype = ET_PENTA6;
		else if (tag == "tet4"  ) etype = ET_TET4;
		else if (tag == "tet10" ) etype = ET_TET10;
		else if (tag == "tet15" ) etype = ET_TET15;
		else if (tag == "quad4" ) etype = ET_QUAD4;
		else if (tag == "tri3"  ) etype = ET_TRI3;
		else if (tag == "truss2") etype = ET_TRUSS2;
		else throw XMLReader::InvalidTag(tag);

		switch (etype)
		{
		case ET_HEX8:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(ne), m_pim->m_nhex8, nid, nmat);
			}
			break;
		case ET_HEX20:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(ne), FE_HEX20G27, nid, nmat);
			}
			break;
		case ET_HEX27:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(ne), FE_HEX27G27, nid, nmat);
			}
			break;
		case ET_PENTA6:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(ne), FE_PENTA6G6, nid, nmat);
			}
			break;
		case ET_TET4:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(ne), m_pim->m_ntet4, nid, nmat);
			}
			break;
		case ET_TET10:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(ne), m_pim->m_ntet10, nid, nmat);
			}
			break;
		case ET_TET15:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(ne), m_pim->m_ntet15, nid, nmat);
			}
			break;
		case ET_QUAD4:
			{
				FEShellDomain& sd = static_cast<FEShellDomain&>(dom);
				ReadShellElement(tag, sd.Element(ne), FE_SHELL_QUAD, nid, nmat);
			}
			break;
		case ET_TRI3:
			{
				FEShellDomain& sd = static_cast<FEShellDomain&>(dom);
				ReadShellElement(tag, sd.Element(ne), FE_SHELL_TRI, nid, nmat);
			}
			break;
		case ET_TRUSS2:
			{
				FETrussDomain& td = static_cast<FETrussDomain&>(dom);
				ReadTrussElement(tag, td.Element(ne), FE_TRUSS, nid, nmat);
			}
			break;
		default:
			throw FEBioImport::InvalidElementType();
		}

		// go to next tag
		++tag;
	}

	// assign material point data
	for (i=0; i<mesh.Domains(); ++i)
	{
		FEDomain& d = mesh.Domain(i);
		d.InitMaterialPointData();
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
	const char* szname = tag.AttributeValue("elset", true);

	// get the element type
	FE_Element_Shape etype;
	const char* sztype = tag.AttributeValue("type");
	if      (strcmp(sztype, "hex8"  ) == 0) etype = ET_HEX8;
	else if (strcmp(sztype, "hex20" ) == 0) etype = ET_HEX20;
	else if (strcmp(sztype, "hex27" ) == 0) etype = ET_HEX27;
	else if (strcmp(sztype, "penta6") == 0) etype = ET_PENTA6;
	else if (strcmp(sztype, "tet4"  ) == 0) etype = ET_TET4;
	else if (strcmp(sztype, "tet10" ) == 0) etype = ET_TET10;
	else if (strcmp(sztype, "tet15" ) == 0) etype = ET_TET15;
	else if (strcmp(sztype, "quad4" ) == 0) etype = ET_QUAD4;
	else if (strcmp(sztype, "tri3"  ) == 0) etype = ET_TRI3;
	else if (strcmp(sztype, "truss2") == 0) etype = ET_TRUSS2;
	else 
	{
		// new way for defining element type and integration rule at the same time
		// this is useful for multi-step analyses where the geometry is read in before the control section.
		if      (strcmp(sztype, "TET10G4"  ) == 0) { etype = ET_TET10; m_pim->m_ntet10 = FE_TET10G4;  }
		else if (strcmp(sztype, "TET10G8"  ) == 0) { etype = ET_TET10; m_pim->m_ntet10 = FE_TET10G8;  }
		else if (strcmp(sztype, "TET10GL11") == 0) { etype = ET_TET10; m_pim->m_ntet10 = FE_TET10GL11;}
		else if (strcmp(sztype, "TET10G4_S3"  ) == 0) { etype = ET_TET10; m_pim->m_ntet10 = FE_TET10G4  ; m_pim->m_ntri6 = FE_TRI6G3; }
		else if (strcmp(sztype, "TET10G8_S3"  ) == 0) { etype = ET_TET10; m_pim->m_ntet10 = FE_TET10G8  ; m_pim->m_ntri6 = FE_TRI6G3; }
		else if (strcmp(sztype, "TET10GL11_S3") == 0) { etype = ET_TET10; m_pim->m_ntet10 = FE_TET10GL11; m_pim->m_ntri6 = FE_TRI6G3; }
		else if (strcmp(sztype, "TET10G4_S4"  ) == 0) { etype = ET_TET10; m_pim->m_ntet10 = FE_TET10G4  ; m_pim->m_ntri6 = FE_TRI6G4; }
		else if (strcmp(sztype, "TET10G8_S4"  ) == 0) { etype = ET_TET10; m_pim->m_ntet10 = FE_TET10G8  ; m_pim->m_ntri6 = FE_TRI6G4; }
		else if (strcmp(sztype, "TET10GL11_S4") == 0) { etype = ET_TET10; m_pim->m_ntet10 = FE_TET10GL11; m_pim->m_ntri6 = FE_TRI6G4; }
		else if (strcmp(sztype, "TET10G4_S7"  ) == 0) { etype = ET_TET10; m_pim->m_ntet10 = FE_TET10G4  ; m_pim->m_ntri6 = FE_TRI6G7; }
		else if (strcmp(sztype, "TET10G8_S7"  ) == 0) { etype = ET_TET10; m_pim->m_ntet10 = FE_TET10G8  ; m_pim->m_ntri6 = FE_TRI6G7; }
		else if (strcmp(sztype, "TET10GL11_S7") == 0) { etype = ET_TET10; m_pim->m_ntet10 = FE_TET10GL11; m_pim->m_ntri6 = FE_TRI6G7; }
		else if (strcmp(sztype, "TET15G8"  ) == 0) { etype = ET_TET15; m_pim->m_ntet15 = FE_TET15G8;  }
		else if (strcmp(sztype, "TET15G11" ) == 0) { etype = ET_TET15; m_pim->m_ntet15 = FE_TET15G11; }
		else if (strcmp(sztype, "TET15G15" ) == 0) { etype = ET_TET15; m_pim->m_ntet15 = FE_TET15G15; }
		else if (strcmp(sztype, "TET15G8_S3"  ) == 0) { etype = ET_TET15; m_pim->m_ntet15 = FE_TET15G8 ; m_pim->m_ntri7 = FE_TRI7G3;}
		else if (strcmp(sztype, "TET15G11_S3" ) == 0) { etype = ET_TET15; m_pim->m_ntet15 = FE_TET15G11; m_pim->m_ntri7 = FE_TRI7G3;}
		else if (strcmp(sztype, "TET15G15_S3" ) == 0) { etype = ET_TET15; m_pim->m_ntet15 = FE_TET15G15; m_pim->m_ntri7 = FE_TRI7G3;}
		else if (strcmp(sztype, "TET15G8_S4"  ) == 0) { etype = ET_TET15; m_pim->m_ntet15 = FE_TET15G8 ; m_pim->m_ntri7 = FE_TRI7G4;}
		else if (strcmp(sztype, "TET15G11_S4" ) == 0) { etype = ET_TET15; m_pim->m_ntet15 = FE_TET15G11; m_pim->m_ntri7 = FE_TRI7G4;}
		else if (strcmp(sztype, "TET15G15_S4" ) == 0) { etype = ET_TET15; m_pim->m_ntet15 = FE_TET15G15; m_pim->m_ntri7 = FE_TRI7G4;}
		else if (strcmp(sztype, "TET15G8_S7"  ) == 0) { etype = ET_TET15; m_pim->m_ntet15 = FE_TET15G8 ; m_pim->m_ntri7 = FE_TRI7G7;}
		else if (strcmp(sztype, "TET15G11_S7" ) == 0) { etype = ET_TET15; m_pim->m_ntet15 = FE_TET15G11; m_pim->m_ntri7 = FE_TRI7G7;}
		else if (strcmp(sztype, "TET15G15_S7" ) == 0) { etype = ET_TET15; m_pim->m_ntet15 = FE_TET15G15; m_pim->m_ntri7 = FE_TRI7G7;}
		else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
	}

	// get the domain's material class
	FEMaterial* pmat = fem.GetMaterial(nmat);

	// then, find the domain type depending on the 
	// element and material types
	int ndomtype = DomainType(etype, pmat);
	if (ndomtype == 0) throw FEBioImport::InvalidDomainType();

	// create the new domain
	FEDomain* pdom = CreateDomain(ndomtype, &mesh, pmat);
	if (pdom == 0) throw FEBioImport::FailedCreatingDomain();
	FEDomain& dom = *pdom;
	dom.SetName(szname);

	// count elements
	int elems = tag.children();
	assert(elems);

	// add domain it to the mesh
	pdom->create(elems);
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
	int nid = m_pim->m_maxid + 1;
	for (int i=0; i<elems; ++i, ++nid)
	{
		if ((tag == "elem")==false) throw XMLReader::InvalidTag(tag);

		// keep track of the largest element ID
		if (nid > m_pim->m_maxid) m_pim->m_maxid = nid;

		// add to the element set (if we have one)
		if (pg) (*pg)[i] = nid;

		// read the element data
		switch (etype)
		{
		case ET_HEX8:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(i), m_pim->m_nhex8, nid, nmat);
			}
			break;
		case ET_PENTA6:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(i), FE_PENTA6G6, nid, nmat);
			}
			break;
		case ET_TET4:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(i), m_pim->m_ntet4, nid, nmat);
			}
			break;
		case ET_TET10:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(i), m_pim->m_ntet10, nid, nmat);
			}
			break;
		case ET_TET15:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(i), m_pim->m_ntet15, nid, nmat);
			}
			break;
		case ET_HEX20:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(i), FE_HEX20G27, nid, nmat);
			}
			break;
		case ET_HEX27:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(i), FE_HEX27G27, nid, nmat);
			}
			break;
		case ET_QUAD4:
			{
				FEShellDomain& sd = static_cast<FEShellDomain&>(dom);
				ReadShellElement(tag, sd.Element(i), FE_SHELL_QUAD, nid, nmat);
			}
			break;
		case ET_TRI3:
			{
				FEShellDomain& sd = static_cast<FEShellDomain&>(dom);
				ReadShellElement(tag, sd.Element(i), FE_SHELL_TRI, nid, nmat);
			}
			break;
		case ET_TRUSS2:
			{
				FETrussDomain& td = static_cast<FETrussDomain&>(dom);
				ReadTrussElement(tag, td.Element(i), FE_TRUSS, nid, nmat);
			}
			break;
		default:
			throw FEBioImport::InvalidElementType();
		}

		// go to next tag
		++tag;
	}

	// assign material point data
	dom.InitMaterialPointData();
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
		FE_Element_Shape etype = ElementShape(t);

		// find a domain for this element
		int ndom = -1;
		for (i=0; i<(int) dom.size(); ++i)
		{
			FEDOMAIN& d = dom[i];
			if (d.elem == etype)
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
			d.elem = etype;
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

		// then, find the domain type depending on the 
		// element and material types
		int ntype = DomainType(d.elem, 0);
		if (ntype == 0) throw FEBioImport::InvalidDomainType();

		// create the new domain
		FEDomain* pdom = CreateDomain(ntype, &mesh, 0);
		if (pdom == 0) throw FEBioImport::FailedCreatingDomain();


		// add it to the mesh
		assert(d.nel);
		pdom->create(d.nel);
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
		int etype = -1;
		if      (tag == "hex8"  ) etype = ET_HEX8;
		else if (tag == "hex20" ) etype = ET_HEX20;
		else if (tag == "hex27" ) etype = ET_HEX27;
		else if (tag == "penta6") etype = ET_PENTA6;
		else if (tag == "tet4"  ) etype = m_pim->m_ntet4;
		else if (tag == "tet10" ) etype = ET_TET10;
		else if (tag == "tet15" ) etype = ET_TET15;
		else if (tag == "quad4" ) etype = ET_QUAD4;
		else if (tag == "tri3"  ) etype = ET_TRI3;
		else if (tag == "truss2") etype = ET_TRUSS2;
		else throw XMLReader::InvalidTag(tag);

		switch (etype)
		{
		case ET_HEX8:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(ne), m_pim->m_nhex8, nid, 0);
			}
			break;
		case ET_HEX20:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(ne), FE_HEX20G27, nid, 0);
			}
			break;
		case ET_HEX27:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(ne), FE_HEX27G27, nid, 0);
			}
			break;
		case ET_PENTA6:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(ne), FE_PENTA6G6, nid, 0);
			}
			break;
		case ET_TET4:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(ne), m_pim->m_ntet4, nid, 0);
			}
			break;
		case ET_TET10:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(ne), FE_TET10G4, nid, 0);
			}
			break;
		case ET_TET15:
			{
				FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
				ReadSolidElement(tag, bd.Element(ne), FE_TET15G8, nid, 0);
			}
			break;
		case ET_QUAD4:
			{
				FEShellDomain& sd = static_cast<FEShellDomain&>(dom);
				ReadShellElement(tag, sd.Element(ne), FE_SHELL_QUAD, nid, 0);
			}
			break;
		case ET_TRI3:
			{
				FEShellDomain& sd = static_cast<FEShellDomain&>(dom);
				ReadShellElement(tag, sd.Element(ne), FE_SHELL_TRI, nid, 0);
			}
			break;
		case ET_TRUSS2:
			{
				FETrussDomain& td = static_cast<FETrussDomain&>(dom);
				ReadTrussElement(tag, td.Element(ne), FE_TRUSS, nid, 0);
			}
			break;
		default:
			throw FEBioImport::InvalidElementType();
		}

		// go to next tag
		++tag;
	}
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection::ReadSolidElement(XMLTag &tag, FESolidElement& el, int ntype, int nid, int nmat)
{
	el.SetType(ntype);
	el.m_nID = nid;
	int n[FEElement::MAX_NODES];
	tag.value(n,el.Nodes());
	for (int j=0; j<el.Nodes(); ++j) el.m_node[j] = n[j]-1;
	el.SetMatID(nmat);
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection::ReadShellElement(XMLTag &tag, FEShellElement& el, int ntype, int nid, int nmat)
{
	el.SetType(ntype);
	el.m_nID = nid;
	int n[9];
	tag.value(n,el.Nodes());
	for (int j=0; j<el.Nodes(); ++j) { el.m_node[j] = n[j]-1; el.m_h0[j] = 0.0; }
	el.SetMatID(nmat);
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection::ReadTrussElement(XMLTag &tag, FETrussElement& el, int ntype, int nid, int nmat)
{
	el.SetType(ntype);
	el.m_nID = nid;
	int n[8];
	tag.value(n, el.Nodes());
	for (int j=0; j<el.Nodes(); ++j) el.m_node[j] = n[j]-1;
	el.SetMatID(nmat);

	// area is read in the ElementData section
	el.m_a0 = 0;
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
			assert(pelem[el.m_nID-1] == 0);
			pelem[el.m_nID-1] = &el;
		}
	}

	// read additional element data
	++tag;
	do
	{
		// make sure this is an "element" tag
		if (tag != "element") throw XMLReader::InvalidTag(tag);

		// get the element number
		const char* szid = tag.AttributeValue("id");
		int n = atoi(szid)-1;

		// make sure the number is valid
		if ((n<0) || (n>=nelems)) throw XMLReader::InvalidAttributeValue(tag, "id", szid);

		// get a pointer to the element
		FEElement* pe = pelem[n];

		vec3d a;
		++tag;
		do
		{
			if (tag == "fiber")
			{
				// read the fiber direction
				m_pim->value(tag, a);

				// normalize fiber
				a.unit();

				// set up a orthonormal coordinate system
				vec3d b(0,1,0);
				if (fabs(fabs(a*b) - 1) < 1e-7) b = vec3d(0,0,1);
				vec3d c = a^b;
				b = c^a;

				// make sure they are unit vectors
				b.unit();
				c.unit();

				for (int i=0; i<pe->GaussPoints(); ++i)
				{
					FEElasticMaterialPoint& pt = *pe->GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
					mat3d& m = pt.m_Q;
					m.zero();
					m[0][0] = a.x; m[0][1] = b.x; m[0][2] = c.x;
					m[1][0] = a.y; m[1][1] = b.y; m[1][2] = c.y;
					m[2][0] = a.z; m[2][1] = b.z; m[2][2] = c.z;
				}
			}
			else if (tag == "mat_axis")
			{
				vec3d a, d;

				++tag;
				do
				{
					if      (tag == "a") m_pim->value(tag, a);
					else if (tag == "d") m_pim->value(tag, d);
					else throw XMLReader::InvalidTag(tag);

					++tag;
				}
				while (!tag.isend());

				vec3d c = a^d;
				vec3d b = c^a;

				// normalize
				a.unit();
				b.unit();
				c.unit();

				// assign to element
				for (int i=0; i<pe->GaussPoints(); ++i)
				{
					FEElasticMaterialPoint& pt = *pe->GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
					mat3d& m = pt.m_Q;
					m.zero();
					m[0][0] = a.x; m[0][1] = b.x; m[0][2] = c.x;
					m[1][0] = a.y; m[1][1] = b.y; m[1][2] = c.y;
					m[2][0] = a.z; m[2][1] = b.z; m[2][2] = c.z;
				}
			}
			else if (tag == "thickness")
			{
				if (pe->Class() != FE_ELEM_SHELL) throw XMLReader::InvalidTag(tag);
				FEShellElement* pse = static_cast<FEShellElement*> (pe);

				// read shell thickness
				tag.value(&pse->m_h0[0],pse->Nodes());
			}
			else if (tag == "area")
			{
				if (pe->Class() != FE_ELEM_TRUSS) throw XMLReader::InvalidTag(tag);
				FETrussElement* pt = static_cast<FETrussElement*>(pe);

				// read truss area
				m_pim->value(tag, pt->m_a0);
			}
			else
			{
				for (int i=0; i<pe->GaussPoints(); ++i)
				{
					FEMaterialPoint* pt = pe->GetMaterialPoint(i);
					while (pt)
					{
						FEParameterList& pl = pt->GetParameterList();
						if (m_pim->ReadParameter(tag, pl)) break;
						else
						{
							pt = pt->Next();
							if (pt == 0) throw XMLReader::InvalidTag(tag);
						}
					}
				}
			}
			++tag;
		}
		while (!tag.isend());

		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
//! Reads the Geometry::Groups section of the FEBio input file

void FEBioGeometrySection::ParseNodeSetSection(XMLTag& tag)
{
	FEMesh& mesh = *m_pim->GetFEMesh();

	// get the name attribute
	const char* szname = tag.AttributeValue("name");

	// create a new node set
	FENodeSet* pns = new FENodeSet(&mesh);
	pns->SetName(szname);

	// NOTE: The initial specs of the 2.0 format defined the nodes as a value list.
	// This was later revised to promote more consistency between the node set definitions
	// and other features (e.g. fixed BC's) that define node sets. However, for now
	// we still support the old format, although this is to be considered obsolete.
	vector<int> l;
	if (tag.isleaf())
	{
		// NOTE: This is the old way of specifying node sets. Consider obsolete!
		// read the list of indices
		m_pim->ReadList(tag, l);
		assert(!l.empty());
	}
	else
	{
		// This is new and preferred way of defining node sets
		++tag;
		do
		{
			if (tag == "node")
			{
				int nid = -1;
				tag.AttributeValue("id", nid);
				l.push_back(nid);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}

	// only add non-empty node sets
	if (l.empty() == false)
	{
		// assign indices to node set
		int N = l.size();
		pns->create(N);
		for (int i=0; i<N; ++i) (*pns)[i] = l[i] - 1;

		// add the nodeset to the mesh
		mesh.AddNodeSet(pns);
	}
	else delete pns;
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
