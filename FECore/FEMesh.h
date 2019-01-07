#pragma once
#include "FENode.h"
#include "FENodeElemList.h"
#include "DumpStream.h"
#include "FENodeSet.h"
#include "FEFacetSet.h"
#include "FEDiscreteSet.h"
#include "FESegmentSet.h"
#include "FEElementSet.h"
#include "FESurfacePair.h"
#include "FEBoundingBox.h"

//-----------------------------------------------------------------------------
class FEEdge;
class FESurface;
class FEDomain;
class FEModel;
class FETimeInfo;

//---------------------------------------------------------------------------------------
// Helper class for faster lookup of elements based on their ID 
class FECORE_API FEElementLUT
{
public:
	FEElementLUT(FEMesh& mesh);

	// Find an element from its ID
	FEElement* Find(int nid);

private:
	vector<FEElement*>	m_elem;
	int					m_minID, m_maxID;
};

//-----------------------------------------------------------------------------
//! Defines a finite element mesh

//! All the geometry data is stored in this class. 

class FECORE_API FEMesh
{
public:
	//! constructor
	FEMesh(FEModel* fem);

	//! destructor
	virtual ~FEMesh();

	//! stream mesh data
	void Serialize(DumpStream& dmp);

    //! initialize material points in mesh
    void InitMaterialPoints();
    
	//! clear the mesh
	void Clear();

	//! allocate storage for mesh data
	void CreateNodes(int nodes);
	void AddNodes(int nodes);

	//! return number of nodes
	int Nodes() const;

	//! return total nr of elements
	int Elements() const;

	//! return the nr of elements of a specific domain type
	int Elements(int ndom_type) const;

	//! return reference to a node
	FENode& Node(int i);
	const FENode& Node(int i) const;

	//! Set the number of degrees of freedom on this mesh
	void SetDOFS(int n);

	//! update bounding box
	void UpdateBox();

	//! retrieve the bounding box
	FEBoundingBox& GetBoundingBox() { return m_box; }

	//! remove isolated vertices
	int RemoveIsolatedVertices();

	//! Reset the mesh data
	void Reset();

	//! Calculates an elements volume
	double ElementVolume(FEElement& el);

	//! Get the face nodes from a given element
	int GetFace(FEElement& el, int n, int* nf);

	//! return the nr of faces an element has
	int Faces(FEElement& el);

	//! Finds a node from a given ID
	FENode* FindNodeFromID(int nid);

	//! return an element (expensive way!)
	FEElement* Element(int i);

	//! Finds an element from a given ID
	FEElement* FindElementFromID(int nid);

	//! Finds the solid element in which y lies
	FESolidElement* FindSolidElement(vec3d y, double r[3]);

	FENodeElemList& NodeElementList()
	{
		if (m_NEL.Size() != m_Node.size()) m_NEL.Create(*this);
		return m_NEL;
	}

	// --- NODESETS ---
	//! adds a node set to the mesh
	void AddNodeSet(FENodeSet* pns) { m_NodeSet.push_back(pns); }

	//! number of nodesets
	int NodeSets() { return (int) m_NodeSet.size(); }

	//! return a node set
	FENodeSet* NodeSet(int i) { return m_NodeSet[i]; }

	//! Find a nodeset by ID
	FENodeSet* FindNodeSet(int nid);

	//! Find a nodeset by name
	FENodeSet* FindNodeSet(const std::string& name);

	// --- ELEMENT SETS ---
	int ElementSets() { return (int) m_ElemSet.size(); }
	FEElementSet& ElementSet(int n) { return *m_ElemSet[n]; }
	void AddElementSet(FEElementSet* pg) { m_ElemSet.push_back(pg); }

	//! Find a element set by name
	FEElementSet* FindElementSet(const std::string& name);

	// --- DOMAINS ---
	int Domains();
	FEDomain& Domain(int n);

	void AddDomain(FEDomain* pd);

	FEDomain* FindDomain(const std::string& name);

	//! get a list of domains that belong to a specific material
	void DomainListFromMaterial(vector<int>& lmat, vector<int>& ldom);

	// --- SURFACES ---
	int Surfaces() { return (int) m_Surf.size(); }
	FESurface& Surface(int n) { return *m_Surf[n]; }
	void AddSurface(FESurface* ps) { m_Surf.push_back(ps); }
	FESurface* FindSurface(const std::string& name);

	// --- EDGES ---
	int Edges() { return (int) m_Edge.size(); }
	FEEdge& Edge(int n) { return *m_Edge[n]; }
	void AddEdge(FEEdge* ps) { m_Edge.push_back(ps); }

	// --- FACETSETS ---
	int FacetSets() { return (int) m_FaceSet.size(); }
	FEFacetSet& FacetSet(int n) { return *m_FaceSet[n]; }
	void AddFacetSet(FEFacetSet* ps) { m_FaceSet.push_back(ps); }
	FEFacetSet* FindFacetSet(const std::string& name);

	// --- Segment Sets ---
	int SegmentSets() { return (int) m_LineSet.size(); }
	FESegmentSet& SegmentSet(int n) { return *m_LineSet[n]; }
	void AddSegmentSet(FESegmentSet* ps) { m_LineSet.push_back(ps); }
	FESegmentSet* FindSegmentSet(const std::string& name);

	// --- Discrete Element Sets ---
	int DiscreteSets() { return (int) m_DiscSet.size(); }
	FEDiscreteSet& DiscreteSet(int n) { return *m_DiscSet[n]; }
	void AddDiscreteSet(FEDiscreteSet* ps) { m_DiscSet.push_back(ps); }
	FEDiscreteSet* FindDiscreteSet(const std::string& name);

	// --- surface pairs ---
	int SurfacePairs() { return (int)m_SurfPair.size(); }
	FESurfacePair& SurfacePair(int n) { return *m_SurfPair[n]; }
	void AddSurfacePair(FESurfacePair* ps) { m_SurfPair.push_back(ps); }
	FESurfacePair* FindSurfacePair(const std::string& name);

public:
	//! Calculate the surface representing the element boundaries
	//! boutside : include all exterior facets
	//! binside  : include all interior facets
	FESurface* ElementBoundarySurface(bool boutside = true, bool binside = false);

	//! Calculate the surface representing the element boundaries
	//! domains  : a list of which domains to create the surface from 
	//! boutside : include all exterior facets
	//! binside  : include all interior facets
	FESurface* ElementBoundarySurface(std::vector<FEDomain*> domains, bool boutside = true, bool binside = false);

	//! get the nodal coordinates in reference configuration
	void GetInitialNodalCoordinates(const FEElement& el, vec3d* node);

	//! get the nodal coordinates in current configuration
	void GetNodalCoordinates(const FEElement& el, vec3d* node);

	// Get the FE model
	FEModel* GetFEModel() const { return m_fem; }

	// update the domains of the mesh
	void Update(const FETimeInfo& tp);

protected:
	double SolidElementVolume(FESolidElement& el);
	double ShellElementVolume(FEShellElement& el);

private:
	vector<FENode>		m_Node;		//!< nodes
	vector<FEDomain*>	m_Domain;	//!< list of domains
	vector<FESurface*>	m_Surf;		//!< surfaces
	vector<FEEdge*>		m_Edge;		//!< Edges

	vector<FENodeSet*>		m_NodeSet;	//!< node sets
	vector<FESegmentSet*>	m_LineSet;	//!< segment sets
	vector<FEFacetSet*>		m_FaceSet;	//!< facet sets
	vector<FEElementSet*>	m_ElemSet;	//!< element sets
	vector<FEDiscreteSet*>	m_DiscSet;	//!< discrete element sets
	vector<FESurfacePair*>	m_SurfPair;	//!< facet set pairs

	FEBoundingBox		m_box;	//!< bounding box

	FENodeElemList	m_NEL;
	FEElementLUT*	m_LUT;

	FEModel*	m_fem;
private:
	//! hide the copy constructor
	FEMesh(FEMesh& m){}

	//! hide assignment operator
	void operator =(FEMesh& m) {}
};
