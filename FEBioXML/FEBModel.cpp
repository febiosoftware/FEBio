#include "stdafx.h"
#include "FEBModel.h"
#include <FECore/FEModel.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FEMaterial.h>

//=============================================================================
FEBModel::NodeSet::NodeSet() {}

FEBModel::NodeSet::NodeSet(const FEBModel::NodeSet& set)
{
	m_name = set.m_name;
	m_node = set.m_node;
}

FEBModel::NodeSet::NodeSet(const string& name) : m_name(name) {}

void FEBModel::NodeSet::SetName(const string& name) { m_name = name; }

const string& FEBModel::NodeSet::Name() const { return m_name; }

void FEBModel::NodeSet::SetNodeList(const vector<int>& node) { m_node = node; }

const vector<int>& FEBModel::NodeSet::NodeList() const { return m_node; }

//=============================================================================
FEBModel::ElementSet::ElementSet() {}

FEBModel::ElementSet::ElementSet(const FEBModel::ElementSet& set)
{
	m_name = set.m_name;
	m_elem = set.m_elem;
}

FEBModel::ElementSet::ElementSet(const string& name) : m_name(name) {}

void FEBModel::ElementSet::SetName(const string& name) { m_name = name; }

const string& FEBModel::ElementSet::Name() const { return m_name; }

void FEBModel::ElementSet::SetElementList(const vector<int>& elem) { m_elem = elem; }

const vector<int>& FEBModel::ElementSet::ElementList() const { return m_elem; }

//=============================================================================
FEBModel::Domain::Domain() {}

FEBModel::Domain::Domain(const FEBModel::Domain& dom)
{
	m_spec = dom.m_spec;
	m_name = dom.m_name;
	m_matName = dom.m_matName;
	m_Elem = dom.m_Elem;
}

FEBModel::Domain::Domain(const FE_Element_Spec& spec) : m_spec(spec) {}

const string& FEBModel::Domain::Name() const { return m_name; }

void FEBModel::Domain::SetName(const string& name) { m_name = name; }

void FEBModel::Domain::SetMaterialName(const string& name) { m_matName = name; }

const string& FEBModel::Domain::MaterialName() const { return m_matName; }

void FEBModel::Domain::SetElementList(const vector<ELEMENT>& el) { m_Elem = el; }

const vector<FEBModel::ELEMENT>& FEBModel::Domain::ElementList() const { return m_Elem; }

//=============================================================================
FEBModel::Surface::Surface() {}

FEBModel::Surface::Surface(const FEBModel::Surface& surf)
{
	m_name = surf.m_name;
	m_Face = surf.m_Face;
}

FEBModel::Surface::Surface(const string& name) : m_name(name) {}

const string& FEBModel::Surface::Name() const { return m_name; }

void FEBModel::Surface::SetName(const string& name) { m_name = name; }

void FEBModel::Surface::SetFacetList(const vector<FEBModel::FACET>& face) { m_Face = face; }

const vector<FEBModel::FACET>& FEBModel::Surface::FacetList() const { return m_Face; }

//=============================================================================
FEBModel::Part::Part() {}

FEBModel::Part::Part(const std::string& name) : m_name(name) {}

FEBModel::Part::Part(const FEBModel::Part& part)
{
	m_name = part.m_name;
	m_Node = part.m_Node;
	for (size_t i=0; i<part.m_Dom.size() ; ++i) AddDomain (new Domain (*part.m_Dom[i]));
	for (size_t i=0; i<part.m_Surf.size(); ++i) AddSurface(new Surface(*part.m_Surf[i]));
	for (size_t i=0; i<part.m_NSet.size(); ++i) AddNodeSet(new NodeSet(*part.m_NSet[i]));
	for (size_t i=0; i<part.m_ESet.size(); ++i) AddElementSet(new ElementSet(*part.m_ESet[i]));
}

FEBModel::Part::~Part()
{
	for (size_t i=0; i<m_NSet.size(); ++i) delete m_NSet[i];
	for (size_t i=0; i<m_Dom.size(); ++i) delete m_Dom[i];
	for (size_t i = 0; i<m_Surf.size(); ++i) delete m_Surf[i];
}

void FEBModel::Part::SetName(const std::string& name) {	m_name = name; }

const string& FEBModel::Part::Name() const { return m_name; }

void FEBModel::Part::AddNodes(const std::vector<NODE>& nodes)
{
	size_t N0 = m_Node.size();
	size_t N = nodes.size();
	if (N > 0)
	{
		m_Node.resize(N0 + N);
		for (int i=0; i<N; ++i)
		{
			m_Node[N0 + i] = nodes[i];
		}
	}
}

FEBModel::Domain* FEBModel::Part::FindDomain(const string& name)
{
	for (size_t i = 0; i<m_Dom.size(); ++i)
	{
		Domain* dom = m_Dom[i];
		if (dom->Name() == name) return dom;
	}
	return 0;
}

void FEBModel::Part::AddDomain(FEBModel::Domain* dom) { m_Dom.push_back(dom); }

void FEBModel::Part::AddSurface(FEBModel::Surface* surf) { m_Surf.push_back(surf); }

//=============================================================================
FEBModel::FEBModel()
{
}

FEBModel::~FEBModel()
{
	for (size_t i=0; i<m_Part.size(); ++i) delete m_Part[i];
	m_Part.clear();
}

size_t FEBModel::Parts() const
{
	return m_Part.size();
}

FEBModel::Part* FEBModel::AddPart(const std::string& name)
{
	Part* part = new Part(name);
	m_Part.push_back(part);
	return part;
}

void FEBModel::AddPart(FEBModel::Part* part)
{
	m_Part.push_back(part);
}

FEBModel::Part* FEBModel::FindPart(const string& name)
{
	for (size_t i=0; i<m_Part.size(); ++i)
	{
		Part* p = m_Part[i];
		if (p->Name() == name) return p;
	}

	return 0;
}

bool FEBModel::BuildPart(FEModel& fem, Part& part, const FETransform& T)
{
	// we'll need the kernel for creating domains
	FECoreKernel& febio = FECoreKernel::GetInstance();

	FEMesh& mesh = fem.GetMesh();

	// build node-index lookup table
	int noff = -1, maxID = 0;
	int N0 = mesh.Nodes();
	int NN = part.Nodes();
	for (int i=0; i<NN; ++i)
	{
		int nid = part.GetNode(i).id;
		if ((noff < 0) || (nid < noff)) noff = nid;
		if (nid > maxID) maxID = nid;
	}
	vector<int> NLT(maxID - noff + 1, -1);
	for (int i=0; i<NN; ++i)
	{
		int nid = part.GetNode(i).id - noff;
		NLT[nid] = i + N0;
	}

	// build element-index lookup table
	int eoff = -1; maxID = 0;
	int E0 = mesh.Elements();
	int NDOM = part.Domains();
	for (int i=0; i<NDOM; ++i)
	{
		const Domain& dom = part.GetDomain(i);
		int NE = dom.Elements();
		for (int j=0; j<NE; ++j)
		{
			int eid = dom.GetElement(j).id;
			if ((eoff < 0) || (eid < eoff)) eoff = eid;
			if (eid > maxID) maxID = eid;
		}
	}
	vector<int> ELT(maxID - eoff + 1, -1);
	int ecount = E0;
	for (int i = 0; i<NDOM; ++i)
	{
		const Domain& dom = part.GetDomain(i);
		int NE = dom.Elements();
		for (int j = 0; j<NE; ++j)
		{
			int eid = dom.GetElement(j).id - eoff;
			ELT[eid] = ecount++;
		}
	}

	// create the nodes
	int nid = N0;
	mesh.AddNodes(NN);
	int n = 0;
	for (int j = 0; j<NN; ++j)
	{
		NODE& partNode = part.GetNode(j);
		FENode& meshNode = mesh.Node(N0 + n++);

		meshNode.SetID(++nid);
		meshNode.m_r0 = T.Transform(partNode.r);
		meshNode.m_rt = meshNode.m_r0;
	}
	assert(n == NN);

	// get the part name
	string partName = part.Name();

	// process domains
	int eid = E0;
	for (int i = 0; i<NDOM; ++i)
	{
		const Domain& partDomain = part.GetDomain(i);

		// element count
		int elems = partDomain.Elements();

		// get the element spect
		FE_Element_Spec spec = partDomain.ElementSpec();

		// get the material
		string matName = partDomain.MaterialName();
		FEMaterial* mat = fem.FindMaterial(matName.c_str());
		if (mat == 0) return false;

		// create the domain
		FEDomain* dom = febio.CreateDomain(spec, &mesh, mat);
		if (dom == 0) return false;

		dom->Create(elems, spec.etype);
		dom->SetMatID(mat->GetID());

		// process element data
		for (int j = 0; j<elems; ++j)
		{
			const ELEMENT& domElement = partDomain.GetElement(j);

			FEElement& el = dom->ElementRef(j);
			el.SetID(++eid);

			int ne = el.Nodes();
			for (int n = 0; n<ne; ++n) el.m_node[n] = NLT[domElement.node[n] - noff];
		}

		// add the domain to the mesh
		mesh.AddDomain(dom);

		// initialize material point data
		dom->CreateMaterialPointData();
	}

	// create node sets
	int NSets = part.NodeSets();
	for (int i = 0; i<NSets; ++i)
	{
		NodeSet* set = part.GetNodeSet(i);

		// create a new node set
		FENodeSet* feset = new FENodeSet(&mesh);

		// add the name
		string name = partName + "." + set->Name();
		feset->SetName(name.c_str());

		// copy indices
		vector<int> nodeList = set->NodeList();
		int nn = nodeList.size();
		for (int j=0; j<nn; ++j) nodeList[j] = NLT[nodeList[j] - noff];
		feset->add(nodeList);

		// add it to the mesh
		mesh.AddNodeSet(feset);
	}

	// create surfaces
	int Surfs = part.Surfaces();
	for (int i=0; i<Surfs; ++i)
	{
		Surface* surf = part.GetSurface(i);
		int faces = surf->Facets();

		// create a new facet set
		FEFacetSet* fset = new FEFacetSet;
		string name = partName + "." + surf->Name();
		fset->SetName(name.c_str());

		// copy data
		fset->Create(faces);
		for (int j=0; j<faces; ++j)
		{
			FACET& srcFacet = surf->GetFacet(j);
			FEFacetSet::FACET& face = fset->Face(j);

			face.ntype = srcFacet.ntype;
			int nf = srcFacet.ntype;	// we assume that the type also identifies the number of nodes
			for (int n=0; n<nf; ++n) face.node[n] = NLT[srcFacet.node[n] - noff];
		}

		// add it to the mesh
		mesh.AddFacetSet(fset);
	}

	// create element sets
/*	int ESets = part.ElementSets();
	for (int i=0; i<ESets; ++i)
	{
		ElementSet& eset = *part.GetElementSet(i);
		vector<int> elist = eset.ElementList();

		int ne = (int) elist.size();
		FEElementSet* feset = new FEElementSet(&mesh);
		feset->create(ne);
	}
*/
	return true;
}
