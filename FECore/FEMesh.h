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
#include "FEBoundingBox.h"

//-----------------------------------------------------------------------------
class FESurface;

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
	// Node status flags
	enum Status {
		EXCLUDE     = 1,	// exclude node from analysis
		SHELL       = 2,	// this node belongs to a shell
		RIGID_CLAMP = 4,	// this node should be clamped to a rigid body (only applies to shell nodes)
	};

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

	//! see if status flags are set
	bool HasFlags(unsigned int flags) const { return ((m_nstate & flags) != 0); }

	//! set the status flags
	void SetFlags(unsigned int flags) { m_nstate = flags; }

	//! get the status falgs
	unsigned int Flags() const { return m_nstate; }

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
    
    vec3d   m_d0;   //!< initial director

public:	// rigid body data
	unsigned int	m_nstate;	//!< node state flags
	int				m_rid;		//!< rigid body number

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
	FENodeSet();
	FENodeSet(FEMesh* pm);
	FENodeSet(const FENodeSet& ns);

	void operator = (const FENodeSet& ns);

	void create(int n);

	void add(int n);

	void add(const vector<int>& ns);

	void add(const FENodeSet& ns);

	int size() const { return (int) m_Node.size(); }

	int& operator [] (int i) { return m_Node[i]; }

	const int& operator [] (int i) const { return m_Node[i]; }

	void SetID(int n) { m_nID = n; }
	int GetID() const { return m_nID; }

	void SetName(const char* sz);
	const char* GetName() const { return m_szname; }

	vector<int>& GetNodeList() { return m_Node; }

	FENode* Node(int i);
	const FENode* Node(int i) const;

	void Serialize(DumpStream& ar);

protected:
	int			m_nID;
	char		m_szname[256];
	vector<int>	m_Node;		//!< list of nodes

protected:
	FEMesh*	m_pmesh;
};

//-----------------------------------------------------------------------------
//! Defines a discrete element set (i.e. node-pairs)
class FEDiscreteSet
{
public:
	struct NodePair
	{
		int	n0, n1;
	};

public:
	FEDiscreteSet(FEMesh* pm);
	void create(int n);
	int size() const { return (int)m_pair.size(); }

	void add(int n0, int n1);

	void SetName(const char* sz);
	const char* GetName() const { return m_szname; }

	const NodePair& Element(int i) const { return m_pair[i]; }	

	void Serialize(DumpStream& ar);

private:
	FEMesh*				m_pmesh;
	vector<NodePair>	m_pair;		//!< list of discrete elements
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
	FEFacetSet(FEMesh* mesh);
	const char* GetName() { return m_szname; }
	void SetName(const char* sz);

	void Create(int n);

	int Faces() const { return (int) m_Face.size(); }
	FACET& Face(int i);
	const FACET& Face(int i) const;

	void Add(FEFacetSet* pf);

	FENodeSet GetNodeSet();

	void Serialize(DumpStream& ar);

	const FEMesh* GetMesh() const { return m_mesh; }

private:
	char			m_szname[256];
	vector<FACET>	m_Face;

private:
	FEMesh*			m_mesh;
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
	FESegmentSet(FEMesh* pm);
	const char* GetName() { return m_szname; }
	void SetName(const char* sz);

	void Create(int n);

	int Segments() { return (int) m_Seg.size(); }
	SEGMENT& Segment(int i);

	void Serialize(DumpStream& ar);

private:
	vector<SEGMENT>	m_Seg;
	char			m_szname[256];
	FEMesh*			m_mesh;
};

//-----------------------------------------------------------------------------
// This class defines a set of elements
class FEElementSet
{
public:
	//! constructor
	FEElementSet(FEMesh* pm);

	void create(int n);

	int size() { return (int)m_Elem.size(); }

	int& operator [] (int i) { return m_Elem[i]; }

	void SetName(const char* sz);
	const char* GetName() { return m_szname; }

	void Serialize(DumpStream& ar);

protected:
	char		m_szname[256];	//!< name of element set
	FEMesh*		m_pmesh;		//!< pointer to parent mesh
	vector<int>	m_Elem;			//!< list of elements
};

//-----------------------------------------------------------------------------
class FESurfacePair
{
public:
	FESurfacePair(FEMesh* pm);

	void SetName(const char* szname);
	const char* GetName() const;

	FEFacetSet* GetMasterSurface();
	void SetMasterSurface(FEFacetSet* pf);

	FEFacetSet* GetSlaveSurface();
	void SetSlaveSurface(FEFacetSet* pf);

	void Serialize(DumpStream& ar);

private:
	char	m_szname[256];
	FEFacetSet*	m_master;
	FEFacetSet*	m_slave;
	FEMesh*		m_mesh;
};

//---------------------------------------------------------------------------------------
// Helper class for faster lookup of elements based on their ID 
class FEElementLUT
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

class FECOREDLL_EXPORT FEMesh
{
public:
	// shell formulations
	// todo: Can I move this elsewhere? I want the mesh to be independent of the shell formulation, but
	//       this was the easiest way to merge the old and new shell formulation.
	enum SHELL_FORMULATION {
		NEW_SHELL,
		OLD_SHELL,
        EAS_SHELL,
        ANS_SHELL
	};

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
	int Nodes() const { return (int)m_Node.size(); }

	//! return total nr of elements
	int Elements() const;

	//! return the nr of elements of a specific domain type
	int Elements(int ndom_type) const;

	//! return reference to a node
	FENode& Node(int i) { return m_Node[i]; }
	const FENode& Node(int i) const { return m_Node[i]; }

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

	//! get the default shell formulation
	SHELL_FORMULATION GetShellFormulation();

	//! set the default shell formulation
	void SetShellFormulation(SHELL_FORMULATION shellType);

	//! Get the face nodes from a given element
	int GetFace(FEElement& el, int n, int* nf);

	//! return the nr of faces an element has
	int Faces(FEElement& el);

	//! Finds a node from a given ID
	FENode* FindNodeFromID(int nid);

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
	FEFacetSet* FindFacetSet(const char* szname);

	// --- Segment Sets ---
	int SegmentSets() { return (int) m_LineSet.size(); }
	FESegmentSet& SegmentSet(int n) { return *m_LineSet[n]; }
	void AddSegmentSet(FESegmentSet* ps) { m_LineSet.push_back(ps); }
	FESegmentSet* FindSegmentSet(const char* szname);

	// --- Discrete Element Sets ---
	int DiscreteSets() { return (int) m_DiscSet.size(); }
	FEDiscreteSet& DiscreteSet(int n) { return *m_DiscSet[n]; }
	void AddDiscreteSet(FEDiscreteSet* ps) { m_DiscSet.push_back(ps); }
	FEDiscreteSet* FindDiscreteSet(const char* szname);

	// --- surface pairs ---
	int SurfacePairs() { return (int)m_SurfPair.size(); }
	FESurfacePair& SurfacePair(int n) { return *m_SurfPair[n]; }
	void AddSurfacePair(FESurfacePair* ps) { m_SurfPair.push_back(ps); }
	FESurfacePair* FindSurfacePair(const char* szname);

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

protected:
	double SolidElementVolume(FESolidElement& el);
	double ShellNewElementVolume(FEShellElement& el);
	double ShellOldElementVolume(FEShellElementOld& el);

	//! Initialize shells
	void InitShells();

protected:
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

	SHELL_FORMULATION	m_defaultShell;

private:
	//! hide the copy constructor
	FEMesh(FEMesh& m){}

	//! hide assignment operator
	void operator =(FEMesh& m) {}
};

#endif // !defined(AFX_FEMESH_H__81ABA97F_AD5F_4F1D_8EE9_95B67EBA448E__INCLUDED_)
