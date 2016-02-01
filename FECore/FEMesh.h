// FEMesh.h: interface for the FEMesh class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEMESH_H__81ABA97F_AD5F_4F1D_8EE9_95B67EBA448E__INCLUDED_)
#define AFX_FEMESH_H__81ABA97F_AD5F_4F1D_8EE9_95B67EBA448E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEDomain.h"
#include "FENodeElemList.h"
#include "DumpStream.h"
#include "FEEdge.h"

//-----------------------------------------------------------------------------
//  This class stores the coordinates of a bounding box
//
struct FE_BOUNDING_BOX
{
	vec3d	r0, r1; // coordinates of opposite corners

	// center of box
	vec3d center() { return (r0+r1)*0.5; }

	// dimensions of box
	double width () { return (r1.x-r0.x); }
	double height() { return (r1.y-r0.y); }
	double depth () { return (r1.z-r0.z); }

	// max dimension
	double radius() 
	{ 
		double w = width();
		double h = height();
		double d = depth();

		if ((w>=d) && (w>=h)) return w;
		if ((h>=w) && (h>=d)) return h;
		
		return d;
	}

	void operator += (vec3d r)
	{
		if (r.x < r0.x) r0.x = r.x;
		if (r.y < r0.y) r0.y = r.y;
		if (r.z < r0.z) r0.z = r.z;
		if (r.x > r1.x) r1.x = r.x;
		if (r.y > r1.y) r1.y = r.y;
		if (r.z > r1.z) r1.z = r.z;
	}

	void inflate(double dx, double dy, double dz)
	{
		r0.x -= dx; r1.x += dx;
		r0.y -= dy; r1.y += dy;
		r0.z -= dz; r1.z += dz;
	}

	// check whether a point is inside or not
	bool IsInside(vec3d r) 
	{ 
		return ((r.x>=r0.x) && (r.y>=r0.y) && (r.z>=r0.z) && (r.x<=r1.x) && (r.y<=r1.y) && (r.z<=r1.z));
	}
};

//-----------------------------------------------------------------------------
//! This class defines a finite element node

//! It stores nodal positions and nodal equations numbers and more.
//!
//! The m_ID array will store the equation number for the corresponding
//! degree of freedom. Its values can be (a) non-negative (0 or higher) which
//! gives the equation number in the linear system of equations, (b) -1 if the
//! dof is fixed, and (c) < -1 if the dof corresponds to a prescribed dof. In
//! that case the corresponding equation number is given by -ID-2.

class FENode
{
public:
	//! default constructor
	FENode();

	//! copy constructor
	FENode(const FENode& n);

	//! assignment operator
	FENode& operator = (const FENode& n);

	//! Set the number of DOFS
	void SetDOFS(int n);

	//! Get the nodal ID
	int GetID() const { return m_nID; }

	//! Set the node ID
	void SetID(int n) { m_nID = n; }

protected:
	int		m_nID;	//!< nodal ID
	
public: // geometry data
	vec3d	m_r0;	//!< initial position
	vec3d	m_rt;	//!< current position

	vec3d	m_at;	//!< nodal acceleration

	vec3d	m_rp;	//!< position of node at previous time step
	vec3d	m_vp;	//!< previous velocity
	vec3d	m_ap;	//!< previous acceleration

	vec3d	m_Fr;	//!< nodal reaction forces

public:	// rigid body data
	int		m_rid;		//!< rigid body number
	bool	m_bshell;	//!< does this node belong to a non-rigid shell element?
	bool	m_bexclude;	//!< exclude this node from the analysis

public:
	double get(int n) const { return m_val[n]; }
	void set(int n, double v) { m_val[n] = v; }
	void inc(int n, double v) { m_val[n] += v; }
	void dec(int n, double v) { m_val[n] -= v; }
	vec3d get_vec3d(int i, int j, int k) const { return vec3d(m_val[i], m_val[j], m_val[k]); }
	void set_vec3d(int i, int j, int k, const vec3d& v) { m_val[i] = v.x; m_val[j] = v.y, m_val[k] = v.z; }

public:
	vector<int>		m_BC;	//!< boundary condition array
	vector<int>		m_ID;	//!< nodal equation numbers
	vector<double>	m_val;	//!< nodal DOF values
};

//-----------------------------------------------------------------------------
// Forward declaration of FEMesh class.
class FEMesh;

//-----------------------------------------------------------------------------
//! Defines a node set
//
class FENodeSet
{
public:
	FENodeSet(FEMesh* pm);

	void create(int n);

	int size() { return m_Node.size(); }

	int& operator [] (int i) { return m_Node[i]; }

	void SetID(int n) { m_nID = n; }
	int GetID() { return m_nID; }

	void SetName(const char* sz);
	const char* GetName() { return m_szname; }

	vector<int>& GetNodeList() { return m_Node; }

protected:
	FEMesh*	m_pmesh;
	vector<int>	m_Node;		//!< list of nodes

	int	m_nID;
	char	m_szname[256];
};

//-----------------------------------------------------------------------------
//! This class defines a set of facets. This can be used in the creation of
//! surfaces.
class FEFacetSet
{
public:
	struct FACET
	{
		int	node[FEElement::MAX_NODES];
		int	ntype;	//	3=tri3, 4=quad4, 6=tri6, 7=tri7, 8=quad8
	};

public:
	FEFacetSet();
	const char* GetName() { return m_szname; }
	void SetName(const char* sz);

	void Create(int n);

	int Faces() { return (int) m_Face.size(); }
	FACET& Face(int i);

public:
	vector<FACET>	m_Face;
	char			m_szname[256];
};

//-----------------------------------------------------------------------------
//! This class defines a set of segments. This can be used in the creation of edges.
class FESegmentSet
{
public:
	struct SEGMENT
	{
		int	node[FEElement::MAX_NODES];
		int	ntype;	//	2=line2
	};

public:
	FESegmentSet();
	const char* GetName() { return m_szname; }
	void SetName(const char* sz);

	void Create(int n);

	int Segments() { return (int) m_Seg.size(); }
	SEGMENT& Segment(int i);

public:
	vector<SEGMENT>	m_Seg;
	char			m_szname[256];
};

//-----------------------------------------------------------------------------
// This class defines a set of elements
class FEElementSet
{
public:
	//! constructor
	FEElementSet(FEMesh* pm);

	void create(int n);

	int size() { return m_Elem.size(); }

	int& operator [] (int i) { return m_Elem[i]; }

	void SetName(const char* sz);
	const char* GetName() { return m_szname; }

protected:
	char		m_szname[256];	//!< name of element set
	FEMesh*		m_pmesh;		//!< pointer to parent mesh
	vector<int>	m_Elem;			//!< list of elements
};

//-----------------------------------------------------------------------------
//! Defines a finite element mesh

//! All the geometry data is stored in this class. 

class FEMesh  
{
public:
	//! constructor
	FEMesh();

	//! destructor
	virtual ~FEMesh();

	//! stream mesh data
	void Serialize(DumpStream& dmp);

	//! initialize mesh
	bool Init();

	//! clear the mesh
	void Clear();

	//! allocate storage for mesh data
	void CreateNodes(int nodes);
	void AddNodes(int nodes);

	//! return number of nodes
	int Nodes() { return m_Node.size(); }

	//! return total nr of elements
	int Elements();

	//! return the nr of elements of a specific domain type
	int Elements(int ndom_type);

	//! return reference to a node
	FENode& Node(int i) { return m_Node[i]; }

	//! Set the number of degrees of freedom on this mesh
	void SetDOFS(int n);

	//! update bounding box
	void UpdateBox();

	//! retrieve the bounding box
	FE_BOUNDING_BOX& GetBoundingBox() { return m_box; }

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
	FENodeSet* FindNodeSet(const char* szname);

	// --- ELEMENT SETS ---
	int ElementSets() { return (int) m_ElemSet.size(); }
	FEElementSet& ElementSet(int n) { return *m_ElemSet[n]; }
	void AddElementSet(FEElementSet* pg) { m_ElemSet.push_back(pg); }

	//! Find a element set by name
	FEElementSet* FindElementSet(const char* szname);

	// --- DOMAINS ---
	int Domains() { return (int) m_Domain.size(); }
	FEDomain& Domain(int n) { return *m_Domain[n]; }

	void AddDomain(FEDomain* pd) { m_Domain.push_back(pd); }

	//! get a list of domains that belong to a specific material
	void DomainListFromMaterial(vector<int>& lmat, vector<int>& ldom);

	// --- SURFACES ---
	int Surfaces() { return (int) m_Surf.size(); }
	FESurface& Surface(int n) { return *m_Surf[n]; }
	void AddSurface(FESurface* ps) { m_Surf.push_back(ps); }

	// --- EDGES ---
	int Edges() { return (int) m_Edge.size(); }
	FEEdge& Edge(int n) { return *m_Edge[n]; }
	void AddEdge(FEEdge* ps) { m_Edge.push_back(ps); }

	// --- FACETSETS ---
	int FacetSets() { return (int) m_FaceSet.size(); }
	FEFacetSet& FacetSet(int n) { return *m_FaceSet[n]; }
	void AddFacetSet(FEFacetSet* ps) { m_FaceSet.push_back(ps); }

	// --- Segment Sets ---
	int SegmentSets() { return (int) m_LineSet.size(); }
	FESegmentSet& SegmentSet(int n) { return *m_LineSet[n]; }
	void AddSegmentSet(FESegmentSet* ps) { m_LineSet.push_back(ps); }
	FESegmentSet* FindSegmentSet(const char* szname);

protected:
	double SolidElementVolume(FESolidElement& el);
	double ShellElementVolume(FEShellElement& el);

	//! Initialize shell normals
	void InitShellNormals();

protected:
	vector<FENode>		m_Node;		//!< nodes
	vector<FEDomain*>	m_Domain;	//!< list of domains
	vector<FESurface*>	m_Surf;		//!< surfaces
	vector<FEEdge*>		m_Edge;		//!< Edges

	vector<FENodeSet*>		m_NodeSet;	//!< node sets
	vector<FESegmentSet*>	m_LineSet;	//!< segment sets
	vector<FEFacetSet*>		m_FaceSet;	//!< facet sets
	vector<FEElementSet*>	m_ElemSet;	//!< element sets

	FE_BOUNDING_BOX		m_box;	//!< bounding box

	FENodeElemList	m_NEL;

private:
	//! hide the copy constructor
	FEMesh(FEMesh& m){}

	//! hide assignment operator
	void operator =(FEMesh& m) {}
};

#endif // !defined(AFX_FEMESH_H__81ABA97F_AD5F_4F1D_8EE9_95B67EBA448E__INCLUDED_)
