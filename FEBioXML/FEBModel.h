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



#pragma once
#include <FECore/vec3d.h>
#include <vector>
#include <string>
#include <FECore/FEElement.h>
#include <FECore/FETransform.h>

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
// This is a helper class for parsing the FEBio input file. It manages the geometry
// of the model as composed of parts. After the Geometry section is read in, call the 
// BuildMesh function to generate the mesh for the model.
class FEBModel
{
public:
	struct NODE
	{
		int		id;	// Nodal ID
		vec3d	r;	// Node position
	};

	struct ELEMENT
	{
		int	id;
		int	node[FEElement::MAX_NODES];
	};

	struct FACET
	{
		int id;
		int node[FEElement::MAX_NODES];
		int ntype;
	};

	struct EDGE
	{
		int id;
		int node[FEElement::MAX_NODES];
		int ntype;
	};

	class Domain
	{
	public:
		Domain();
		Domain(const Domain& dom);
		Domain(const FE_Element_Spec& spec);

		void SetName(const std::string& name);
		const std::string& Name() const;

		void SetMaterialName(const std::string& name);
		const std::string& MaterialName() const;

		void SetElementList(const std::vector<ELEMENT>& el);
		const std::vector<ELEMENT>& ElementList() const;

		int Elements() const { return (int) m_Elem.size(); }

		FE_Element_Spec ElementSpec() const { return m_spec; }
		void SetElementSpec(FE_Element_Spec spec) { m_spec = spec; }
		
		void Create(int nsize) { m_Elem.resize(nsize); }
		void Reserve(int nsize) { m_Elem.reserve(nsize); }
		const ELEMENT& GetElement(int i) const { return m_Elem[i]; }
		ELEMENT& GetElement(int i) { return m_Elem[i]; }
		void AddElement(const ELEMENT& el) { m_Elem.push_back(el); }

	private:
		FE_Element_Spec		m_spec;
		std::string			m_name;
		std::string			m_matName;
		std::vector<ELEMENT>	m_Elem;

	public:
		double	m_defaultShellThickness;
	};

	class Surface
	{
	public:
		Surface();
		Surface(const Surface& surf);
		Surface(const std::string& name);

		void SetName(const std::string& name);
		const std::string& Name() const;

		void SetFacetList(const std::vector<FACET>& el);
		const std::vector<FACET>& FacetList() const;

		void Create(int n) { m_Face.resize(n); }

		int Facets() const { return (int) m_Face.size(); }
		FACET& GetFacet(int i) { return m_Face[i]; }

	private:
		std::string	m_name;
		std::vector<FACET>	m_Face;
	};

	class NodeSet
	{
	public:
		NodeSet();
		NodeSet(const NodeSet& set);
		NodeSet(const std::string& name);

		void SetName(const std::string& name);
		const std::string& Name() const;

		void SetNodeList(const std::vector<int>& node);
		const std::vector<int>& NodeList() const;

	private:
		std::string		m_name;
		std::vector<int>	m_node;
	};

	class EdgeSet
	{
	public:
		EdgeSet();
		EdgeSet(const EdgeSet& set);
		EdgeSet(const std::string& name);

		void SetName(const std::string& name);
		const std::string& Name() const;

		void SetEdgeList(const std::vector<EDGE>& edge);
		const std::vector<EDGE>& EdgeList() const;

		int Edges() const { return (int)m_edge.size(); }
		EDGE& Edge(int i) { return m_edge[i]; }

	private:
		std::string			m_name;
		std::vector<EDGE>	m_edge;
	};

	class ElementSet
	{
	public:
		ElementSet();
		ElementSet(const ElementSet& set);
		ElementSet(const std::string& name);

		void SetName(const std::string& name);
		const std::string& Name() const;

		void SetElementList(const std::vector<int>& elem);
		const std::vector<int>& ElementList() const;

	private:
		std::string			m_name;
		std::vector<int>	m_elem;
	};

	class SurfacePair
	{
	public:
		SurfacePair();
		SurfacePair(const SurfacePair& surfPair);

		const std::string& Name() const;

	public:
		std::string	m_name;
		std::string	m_primary;
		std::string	m_secondary;		
	};

	class DiscreteSet
	{
	public:
		struct ELEM
		{
			int	node[2];
		};

	public:
		DiscreteSet();
		DiscreteSet(const DiscreteSet& set);

		void SetName(const std::string& name);
		const std::string& Name() const;

		void AddElement(int n0, int n1);
		const std::vector<ELEM>& ElementList() const;

	private:
		std::string			m_name;
		std::vector<ELEM>	m_elem;
	};

	class Part
	{
	public:
		Part();
		Part(const std::string& name);
		Part(const Part& part);
		~Part();

		void SetName(const std::string& name);
		const std::string& Name() const;

		void AddNodes(const std::vector<NODE>& nodes);

		int Domains() const { return (int)m_Dom.size(); }
		void AddDomain(Domain* dom);
		const Domain& GetDomain(int i) const { return *m_Dom[i]; }
		Domain* FindDomain(const std::string& name);

		int Surfaces() const { return (int) m_Surf.size(); }
		void AddSurface(Surface* surf);
		Surface* GetSurface(int i) { return m_Surf[i]; }
		Surface* FindSurface(const std::string& name);

		int NodeSets() const { return (int) m_NSet.size(); }
		void AddNodeSet(NodeSet* nset) { m_NSet.push_back(nset); }
		NodeSet* GetNodeSet(int i) { return m_NSet[i]; }
		NodeSet* FindNodeSet(const std::string& name);

		int EdgeSets() const { return (int)m_LSet.size(); }
		void AddEdgeSet(EdgeSet* cset) { m_LSet.push_back(cset); }
		EdgeSet* GetEdgeSet(int i) { return m_LSet[i]; }
		EdgeSet* FindEdgeSet(const std::string& name);

		int ElementSets() const { return (int) m_ESet.size(); }
		void AddElementSet(ElementSet* eset) { m_ESet.push_back(eset); }
		ElementSet* GetElementSet(int i) { return m_ESet[i]; }
		ElementSet* FindElementSet(const std::string& name);

		int SurfacePairs() const { return (int)m_SurfPair.size(); }
		void AddSurfacePair(SurfacePair* sp) { m_SurfPair.push_back(sp); }
		SurfacePair* GetSurfacePair(int i) { return m_SurfPair[i]; }

		int DiscreteSets() const { return (int)m_DiscSet.size(); }
		void AddDiscreteSet(DiscreteSet* sp) { m_DiscSet.push_back(sp); }
		DiscreteSet* GetDiscreteSet(int i) { return m_DiscSet[i]; }

		int Nodes() const { return (int) m_Node.size(); }

		NODE& GetNode(int i) { return m_Node[i]; }

	private:
		std::string					m_name;
		std::vector<NODE>			m_Node;
		std::vector<Domain*>		m_Dom;
		std::vector<Surface*>		m_Surf;
		std::vector<NodeSet*>		m_NSet;
		std::vector<EdgeSet*>		m_LSet;
		std::vector<ElementSet*>	m_ESet;
		std::vector<SurfacePair*>	m_SurfPair;
		std::vector<DiscreteSet*>	m_DiscSet;
	};

public:
	FEBModel();
	~FEBModel();

	size_t Parts() const;
	Part* GetPart(int i);
	Part* AddPart(const std::string& name);
	void AddPart(Part* part);

	Part* FindPart(const std::string& name);

	bool BuildPart(FEModel& fem, Part& part, bool buildDomains = true, const Transform& T = Transform());

private:
	std::vector<Part*>	m_Part;
};
