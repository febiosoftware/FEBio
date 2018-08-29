#pragma once
#include <FECore/vec3d.h>
#include <vector>
#include <string>
#include <FECore/FEElement.h>
#include <FECore/FETransform.h>
using namespace std;

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

	class Domain
	{
	public:
		Domain();
		Domain(const Domain& dom);
		Domain(const FE_Element_Spec& spec);

		void SetName(const string& name);
		const string& Name() const;

		void SetMaterialName(const string& name);
		const string& MaterialName() const;

		void SetElementList(const vector<ELEMENT>& el);
		const vector<ELEMENT>& ElementList() const;

		int Elements() const { return (int) m_Elem.size(); }

		FE_Element_Spec ElementSpec() const { return m_spec; }
		
		void Create(int nsize) { m_Elem.resize(nsize); }
		const ELEMENT& GetElement(int i) const { return m_Elem[i]; }
		ELEMENT& GetElement(int i) { return m_Elem[i]; }

	private:
		FE_Element_Spec		m_spec;
		string				m_name;
		string				m_matName;
		vector<ELEMENT>		m_Elem;
	};

	class Surface
	{
	public:
		Surface();
		Surface(const Surface& surf);
		Surface(const string& name);

		void SetName(const string& name);
		const string& Name() const;

		void SetFacetList(const vector<FACET>& el);
		const vector<FACET>& FacetList() const;

		void Create(int n) { m_Face.resize(n); }

		int Facets() const { return (int) m_Face.size(); }
		FACET& GetFacet(int i) { return m_Face[i]; }

	private:
		string	m_name;
		vector<FACET>	m_Face;
	};

	class NodeSet
	{
	public:
		NodeSet();
		NodeSet(const NodeSet& set);
		NodeSet(const string& name);

		void SetName(const string& name);
		const string& Name() const;

		void SetNodeList(const vector<int>& node);
		const vector<int>& NodeList() const;

	private:
		string		m_name;
		vector<int>	m_node;
	};

	class ElementSet
	{
	public:
		ElementSet();
		ElementSet(const ElementSet& set);
		ElementSet(const string& name);

		void SetName(const string& name);
		const string& Name() const;

		void SetElementList(const vector<int>& elem);
		const vector<int>& ElementList() const;

	private:
		string		m_name;
		vector<int>	m_elem;
	};

	class Part
	{
	public:
		Part();
		Part(const std::string& name);
		Part(const Part& part);
		~Part();

		void SetName(const std::string& name);
		const string& Name() const;

		void AddNodes(const std::vector<NODE>& nodes);

		int Domains() const { return (int)m_Dom.size(); }
		void AddDomain(Domain* dom);
		const Domain& GetDomain(int i) const { return *m_Dom[i]; }
		Domain* FindDomain(const string& name);

		int Surfaces() const { return (int) m_Surf.size(); }
		void AddSurface(Surface* surf);
		Surface* GetSurface(int i) { return m_Surf[i]; }

		int NodeSets() const { return (int) m_NSet.size(); }
		void AddNodeSet(NodeSet* nset) { m_NSet.push_back(nset); }
		NodeSet* GetNodeSet(int i) { return m_NSet[i]; }

		int ElementSets() const { return (int) m_ESet.size(); }
		void AddElementSet(ElementSet* eset) { m_ESet.push_back(eset); }
		ElementSet* GetElementSet(int i) { return m_ESet[i]; }

		int Nodes() const { return (int) m_Node.size(); }

		NODE& GetNode(int i) { return m_Node[i]; }

	private:
		string				m_name;
		vector<NODE>		m_Node;
		vector<Domain*>		m_Dom;
		vector<Surface*>	m_Surf;
		vector<NodeSet*>	m_NSet;
		vector<ElementSet*>	m_ESet;
	};

public:
	FEBModel();
	~FEBModel();

	size_t Parts() const;
	Part* AddPart(const std::string& name);
	void AddPart(Part* part);

	Part* FindPart(const string& name);

	bool BuildPart(FEModel& fem, Part& part, const FETransform& T = FETransform());

private:
	std::vector<Part*>	m_Part;
};
