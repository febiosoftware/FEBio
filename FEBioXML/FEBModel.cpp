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
#include "FEBModel.h"
#include <FECore/FEModel.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FEMaterial.h>
#include <FECore/FEDomain.h>
#include <FECore/FEShellDomain.h>
#include <FECore/log.h>
using namespace std;

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
FEBModel::EdgeSet::EdgeSet() {}

FEBModel::EdgeSet::EdgeSet(const FEBModel::EdgeSet& set)
{
	m_name = set.m_name;
	m_edge = set.m_edge;
}

FEBModel::EdgeSet::EdgeSet(const string& name) : m_name(name) {}

void FEBModel::EdgeSet::SetName(const string& name) { m_name = name; }

const string& FEBModel::EdgeSet::Name() const { return m_name; }

void FEBModel::EdgeSet::SetEdgeList(const vector<EDGE>& edge) { m_edge = edge; }

const vector<FEBModel::EDGE>& FEBModel::EdgeSet::EdgeList() const { return m_edge; }

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
FEBModel::SurfacePair::SurfacePair() {}
FEBModel::SurfacePair::SurfacePair(const SurfacePair& surfPair)
{
	m_name = surfPair.m_name;
	m_primary = surfPair.m_primary;
	m_secondary = surfPair.m_secondary;
}
const string& FEBModel::SurfacePair::Name() const { return m_name; }

//=============================================================================
FEBModel::Domain::Domain()
{
	m_defaultShellThickness = 0.0;
}

FEBModel::Domain::Domain(const FEBModel::Domain& dom)
{
	m_spec = dom.m_spec;
	m_name = dom.m_name;
	m_matName = dom.m_matName;
	m_Elem = dom.m_Elem;
	m_defaultShellThickness = dom.m_defaultShellThickness;
}

FEBModel::Domain::Domain(const FE_Element_Spec& spec) : m_spec(spec) 
{
	m_defaultShellThickness = 0.0;
}

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
FEBModel::DiscreteSet::DiscreteSet() {}
FEBModel::DiscreteSet::DiscreteSet(const FEBModel::DiscreteSet& set)
{
	m_name = set.m_name;
	m_elem = set.m_elem;
}

void FEBModel::DiscreteSet::SetName(const string& name) { m_name = name; }
const string& FEBModel::DiscreteSet::Name() const { return m_name; }

void FEBModel::DiscreteSet::AddElement(int n0, int n1) { m_elem.push_back(ELEM{ n0, n1 } ); }
const vector<FEBModel::DiscreteSet::ELEM>& FEBModel::DiscreteSet::ElementList() const { return m_elem; }

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
	for (size_t i = 0; i < part.m_SurfPair.size(); ++i) AddSurfacePair(new SurfacePair(*part.m_SurfPair[i]));
	for (size_t i = 0; i < part.m_DiscSet.size(); ++i) AddDiscreteSet(new DiscreteSet(*part.m_DiscSet[i]));
}

FEBModel::Part::~Part()
{
	for (size_t i=0; i<m_NSet.size(); ++i) delete m_NSet[i];
	for (size_t i=0; i<m_Dom.size(); ++i) delete m_Dom[i];
	for (size_t i = 0; i<m_Surf.size(); ++i) delete m_Surf[i];
	for (size_t i = 0; i < m_SurfPair.size(); ++i) delete m_SurfPair[i];
	for (size_t i = 0; i < m_DiscSet.size(); ++i) delete m_DiscSet[i];
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

FEBModel::Surface* FEBModel::Part::FindSurface(const string& name)
{
	for (size_t i = 0; i < m_Surf.size(); ++i)
	{
		Surface* surf = m_Surf[i];
		if (surf->Name() == name) return surf;
	}
	return nullptr;
}

FEBModel::NodeSet* FEBModel::Part::FindNodeSet(const string& name)
{
	for (size_t i = 0; i < m_NSet.size(); ++i)
	{
		NodeSet* nset = m_NSet[i];
		if (nset->Name() == name) return nset;
	}
	return nullptr;
}

FEBModel::EdgeSet* FEBModel::Part::FindEdgeSet(const string& name)
{
	for (size_t i = 0; i < m_LSet.size(); ++i)
	{
		EdgeSet* lset = m_LSet[i];
		if (lset->Name() == name) return lset;
	}
	return nullptr;
}

FEBModel::ElementSet* FEBModel::Part::FindElementSet(const string& name)
{
	for (size_t i = 0; i < m_ESet.size(); ++i)
	{
		ElementSet* eset = m_ESet[i];
		if (eset->Name() == name) return eset;
	}
	return nullptr;
}

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

FEBModel::Part* FEBModel::GetPart(int i)
{
	return m_Part[i];
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

bool FEBModel::BuildPart(FEModel& fem, Part& part, bool buildDomains, const Transform& T)
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
		meshNode.m_r0 = T.Apply(partNode.r);
		meshNode.m_rt = meshNode.m_r0;
	}
	assert(n == NN);

	// get the part name
	string partName = part.Name();
	if (partName.empty() == false) partName += ".";

	// process domains
	if (buildDomains)
	{
		int eid = E0;
		for (int i = 0; i < NDOM; ++i)
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

			if (dom->Create(elems, spec) == false)
			{
				return false;
			}

			dom->SetMatID(mat->GetID() - 1);

			string domName = partName + partDomain.Name();
			dom->SetName(domName);

			// process element data
			for (int j = 0; j < elems; ++j)
			{
				const ELEMENT& domElement = partDomain.GetElement(j);

				FEElement& el = dom->ElementRef(j);
				el.SetID(++eid);

				int ne = el.Nodes();
				for (int n = 0; n < ne; ++n) el.m_node[n] = NLT[domElement.node[n] - noff];
			}

			if (partDomain.m_defaultShellThickness != 0.0)
			{
				double h0 = partDomain.m_defaultShellThickness;
				FEShellDomain* shellDomain = dynamic_cast<FEShellDomain*>(dom);
				if (shellDomain)
				{
					int ne = shellDomain->Elements();
					for (int n = 0; n < ne; ++n)
					{
						FEShellElement& el = shellDomain->Element(n);
						for (int k = 0; k < el.Nodes(); ++k) el.m_h0[k] = h0;
					}
				}
				else
				{
					FEModel* pfem = &fem;
					feLogWarningEx(pfem, "Shell thickness assigned on non-shell part %s", partDomain.Name().c_str());
				}
			}

			// add the domain to the mesh
			mesh.AddDomain(dom);

			// initialize material point data
			dom->CreateMaterialPointData();
		}
	}

	// create node sets
	int NSets = part.NodeSets();
	for (int i = 0; i<NSets; ++i)
	{
		NodeSet* set = part.GetNodeSet(i);

		// create a new node set
		FENodeSet* feset = new FENodeSet(&fem);

		// add the name
		string name = partName + set->Name();
		feset->SetName(name.c_str());

		// copy indices
		vector<int> nodeList = set->NodeList();
		int nn = (int)nodeList.size();
		for (int j=0; j<nn; ++j) nodeList[j] = NLT[nodeList[j] - noff];
		feset->Add(nodeList);

		// add it to the mesh
		mesh.AddNodeSet(feset);
	}

	// create edges
	int Edges = part.EdgeSets();
	for (int i = 0; i < Edges; ++i)
	{
		EdgeSet* edgeSet = part.GetEdgeSet(i);
		int N = edgeSet->Edges();

		// create a new segment set
		FESegmentSet* segSet = new FESegmentSet(&fem);
		string name = partName + edgeSet->Name();
		segSet->SetName(name.c_str());

		// copy data
		segSet->Create(N);
		for (int j = 0; j < N; ++j)
		{
			EDGE& edge = edgeSet->Edge(j);
			FESegmentSet::SEGMENT& seg = segSet->Segment(j);

			seg.ntype = edge.ntype;
			int nn = edge.ntype;	// we assume that the type also identifies the number of nodes
			for (int n = 0; n < nn; ++n) seg.node[n] = NLT[edge.node[n] - noff];
		}

		// add it to the mesh
		mesh.AddSegmentSet(segSet);
	}

	// create surfaces
	int Surfs = part.Surfaces();
	for (int i=0; i<Surfs; ++i)
	{
		Surface* surf = part.GetSurface(i);
		int faces = surf->Facets();

		// create a new facet set
		FEFacetSet* fset = new FEFacetSet(&fem);
		string name = partName + surf->Name();
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
	int ESets = part.ElementSets();
	for (int i=0; i<ESets; ++i)
	{
		ElementSet& eset = *part.GetElementSet(i);
		vector<int> elist = eset.ElementList();

		int ne = (int) elist.size();
		FEElementSet* feset = new FEElementSet(&fem);
		string name = partName + eset.Name();
		feset->SetName(name);

		// If a domain exists with the same name, we assume
		// that this element set refers to the that domain (TODO: should actually check this!)
		FEDomain* dom = mesh.FindDomain(name);
		if (dom) feset->Create(dom);
		else
		{
			// A domain with the same name is not found, but it is possible that this 
			// set still coincides with a domain, so let's see if we can find it. 
			// see if all elements belong to the same domain
			bool oneDomain = true;
			FEElement* el = mesh.FindElementFromID(elist[0]); assert(el);
			FEDomain* dom = dynamic_cast<FEDomain*>(el->GetMeshPartition());
			for (int i = 1; i < elist.size(); ++i)
			{
				FEElement* el_i = mesh.FindElementFromID(elist[i]); assert(el);
				FEDomain* dom_i = dynamic_cast<FEDomain*>(el_i->GetMeshPartition());

				if (dom != dom_i)
				{
					oneDomain = false;
					break;
				}
			}

			// assign indices to element set
			if (oneDomain && (dom->Elements() == elist.size()))
				feset->Create(dom, elist);
			else
			{
				// Couldn't find a single domain.
				// But maybe this set encompasses the entire mesh? 
				if (elist.size() == mesh.Elements())
				{
					FEDomainList allDomains;
					for (int i = 0; i < mesh.Domains(); ++i) allDomains.AddDomain(&mesh.Domain(i));
					feset->Create(allDomains);
				}
				else
				{
					feset->Create(elist);
				}
			}
		}

		mesh.AddElementSet(feset);
	}

	// create surface pairs
	int SPairs = part.SurfacePairs();
	for (int i = 0; i < SPairs; ++i)
	{
		SurfacePair& spair = *part.GetSurfacePair(i);
		string name = partName + spair.Name();

		FESurfacePair* fesurfPair = new FESurfacePair(&mesh);
		mesh.AddSurfacePair(fesurfPair);
		fesurfPair->SetName(name);

		FEFacetSet* surf1 = mesh.FindFacetSet(spair.m_primary);
		if (surf1 == nullptr) return false;
		fesurfPair->SetPrimarySurface(surf1);

		FEFacetSet* surf2 = mesh.FindFacetSet(spair.m_secondary);
		if (surf2 == nullptr) return false;
		fesurfPair->SetSecondarySurface(surf2);
	}

	// create discrete element sets
	int DSets = part.DiscreteSets();
	for (int i = 0; i < DSets; ++i)
	{
		DiscreteSet& dset = *part.GetDiscreteSet(i);
		string name = partName + dset.Name();

		FEDiscreteSet* fedset = new FEDiscreteSet(&mesh);
		mesh.AddDiscreteSet(fedset);
		fedset->SetName(name);

		const std::vector<DiscreteSet::ELEM>& elemList = dset.ElementList();

		for (int j = 0; j < elemList.size(); ++j)
		{
			int n0 = NLT[elemList[j].node[0] - noff];
			int n1 = NLT[elemList[j].node[1] - noff];

			fedset->add(n0, n1);
		}
	}

	return true;
}
