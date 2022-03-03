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
#include "FENode.h"
#include "FENodeElemList.h"
#include "FENodeSet.h"
#include "FEFacetSet.h"
#include "FEDiscreteSet.h"
#include "FESegmentSet.h"
#include "FEElementSet.h"
#include "FESurfacePair.h"
#include "FEBoundingBox.h"
#include "FESolidElement.h"
#include "FEShellElement.h"

//-----------------------------------------------------------------------------
class FEEdge;
class FESurface;
class FEDomain;
class FEModel;
class FETimeInfo;
class FEDataMap;
class DumpStream;

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

	//! Calculates an elements volume in reference configuration
	double ElementVolume(FEElement& el);

	//! calculates element volume in current configuration
	double CurrentElementVolume(FEElement& el);

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

	//! See if all elements are of a particular shape
	bool IsType(FE_Element_Shape eshape);

	// --- NODESETS ---
	//! adds a node set to the mesh
	void AddNodeSet(FENodeSet* pns) { m_NodeSet.push_back(pns); }

	//! number of nodesets
	int NodeSets() { return (int) m_NodeSet.size(); }

	//! return a node set
	FENodeSet* NodeSet(int i) { return m_NodeSet[i]; }

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
	int FindDomainIndex(const std::string& name);
	FEDomain* FindDomain(int domId);

	//! clear all domains
	void ClearDomains();

	//! Rebuild the LUT
	void RebuildLUT();

	// --- SURFACES ---
	int Surfaces() { return (int) m_Surf.size(); }
	FESurface& Surface(int n) { return *m_Surf[n]; }
	void AddSurface(FESurface* ps) { m_Surf.push_back(ps); }
	FESurface* FindSurface(const std::string& name);
	int FindSurfaceIndex(const std::string& name);

	// create a surface from a facet set
	FESurface* CreateSurface(FEFacetSet& facetSet);

	// --- EDGES ---
	int Edges() { return (int) m_Edge.size(); }
	FEEdge& Edge(int n) { return *m_Edge[n]; }
	void AddEdge(FEEdge* ps) { m_Edge.push_back(ps); }

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

	// --- FACETSETS ---
	int FacetSets() { return (int)m_FaceSet.size(); }
	FEFacetSet& FacetSet(int n) { return *m_FaceSet[n]; }
	void AddFacetSet(FEFacetSet* ps) { m_FaceSet.push_back(ps); }
	FEFacetSet* FindFacetSet(const std::string& name);

	// --- surface pairs ---
	int SurfacePairs() { return (int)m_SurfPair.size(); }
	FESurfacePair& SurfacePair(int n) { return *m_SurfPair[n]; }
	void AddSurfacePair(FESurfacePair* ps) { m_SurfPair.push_back(ps); }
	FESurfacePair* FindSurfacePair(const std::string& name);

public:
	//! stream mesh data
	void Serialize(DumpStream& dmp);

	static void SaveClass(DumpStream& ar, FEMesh* p);
	static FEMesh* LoadClass(DumpStream& ar, FEMesh* p);

	// create a copy of this mesh
	void CopyFrom(FEMesh& mesh);

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

public: // data maps
	void ClearDataMaps();
	void AddDataMap(FEDataMap* map);
	FEDataMap* FindDataMap(const std::string& map);

	int DataMaps() const;
	FEDataMap* GetDataMap(int i);

private:
	vector<FENode>		m_Node;		//!< nodes
	vector<FEDomain*>	m_Domain;	//!< list of domains
	vector<FESurface*>	m_Surf;		//!< surfaces
	vector<FEEdge*>		m_Edge;		//!< Edges

	vector<FENodeSet*>		m_NodeSet;	//!< node sets
	vector<FESegmentSet*>	m_LineSet;	//!< segment sets
	vector<FEElementSet*>	m_ElemSet;	//!< element sets
	vector<FEDiscreteSet*>	m_DiscSet;	//!< discrete element sets
	vector<FEFacetSet*>		m_FaceSet;	//!< facet sets
	vector<FESurfacePair*>	m_SurfPair;	//!< facet set pairs

	vector<FEDataMap*>		m_DataMap;	//!< all data maps

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

class FECORE_API FEElementIterator
{
public:
	FEElementIterator(FEMesh* mesh, FEElementSet* elemSet = nullptr);

	FEElement& operator *() { return *m_el; }

	bool isValid() { return (m_el != nullptr); }

	void operator++();

	void reset();

private:
	FEElement*	m_el;
	FEMesh*	m_mesh;
	FEElementSet*	m_eset;
	int	m_index;
	int	m_dom;
};
