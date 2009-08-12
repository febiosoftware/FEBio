// FESurface.h: interface for the FESurface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FESURFACE_H__6437C4B1_5BB7_4DDA_8354_CADFF3291D3E__INCLUDED_)
#define AFX_FESURFACE_H__6437C4B1_5BB7_4DDA_8354_CADFF3291D3E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEMesh.h"
#include "FENodeElemList.h"
#include "LoadCurve.h"
#include "FENNQuery.h"
#include "mat2d.h"

class FEM;

//-----------------------------------------------------------------------------
//! Surface mesh

//! This class implements the basic functionality for an FE surface.
//! More specialized surfaces are derived from this class

class FESurface  
{
public:
	//! constructor
	FESurface(FEMesh* pmesh) : m_pmesh(pmesh){}

	//! destructor
	virtual ~FESurface(){}

	//! creates surface
	void Create(int n) { m_el.create(n); }

	//! return an element of the surface
	FESurfaceElement& Element(int i) { return m_el[i]; }

	//! return number of surface elements
	int Elements() { return m_el.size(); }

	//! Project a node onto a surface element
	vec3d ProjectToSurface(FESurfaceElement& el, vec3d x, double& r, double& s);

	//! return the mesh to which this surface is attached
	FEMesh* GetMesh() { return m_pmesh; }

	//! number of nodes on this surface
	int Nodes() { return node.size(); }

	//! initialize surface data structure
	virtual void Init();

	//! return the FENode object for local node n
	FENode& Node(int n) { return m_pmesh->Node( node[n] ); }

	//! calculate the surface area of a surface element
	double FaceArea(FESurfaceElement& el);

	//! calculate the metric tensor in the reference configuration
	mat2d Metric0(FESurfaceElement& el, double r, double s);

	//! calculate the surface normal
	vec3d SurfaceNormal(FESurfaceElement& el, double r, double s);

	//! calculate the surface normal at an integration point
	vec3d SurfaceNormal(FESurfaceElement& el, int n);

	//! calculate the global position of a point on the surface
	vec3d Local2Global(FESurfaceElement& el, double r, double s);

	//! calculate the global position of an integration point
	vec3d Local2Global(FESurfaceElement& el, int n);

	//! check to see if a point is on element
	bool IsInsideElement(FESurfaceElement& el, double r, double s, double tol = 0);

	//! return the point on a surface element
	vec3d PointOnSurface(FESurfaceElement& el, double r, double s);

	//! calculates covariant base vectors of a surface
	void CoBaseVectors0(FESurfaceElement& el, double r, double s, vec3d t[2]);

	//! calculates contravariant base vectors of a surface
	void ContraBaseVectors0(FESurfaceElement& el, double r, double s, vec3d t[2]);

protected:
	FEMesh*	m_pmesh;			//!< pointer to parent mesh

	vector<FESurfaceElement>	m_el;	//!< surface elements

public:
	vector<int>	node;	//!< array of node indices

	FENodeElemList	m_NEL;	//!< the node element list
};

#endif // !defined(AFX_FESURFACE_H__6437C4B1_5BB7_4DDA_8354_CADFF3291D3E__INCLUDED_)
