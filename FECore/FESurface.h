// FESurface.h: interface for the FESurface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FESURFACE_H__6437C4B1_5BB7_4DDA_8354_CADFF3291D3E__INCLUDED_)
#define AFX_FESURFACE_H__6437C4B1_5BB7_4DDA_8354_CADFF3291D3E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEMesh.h"
#include "mat2d.h"
#include "vec2d.h"

//-----------------------------------------------------------------------------
//! Surface mesh

//! This class implements the basic functionality for an FE surface.
//! More specialized surfaces are derived from this class

class FESurface : public FEDomain
{
public:
	//! constructor
	FESurface(FEMesh* pm);

	//! destructor
	virtual ~FESurface(){}

	//! initialize surface data structure
	virtual bool Init();

	//! creates surface
	void create(int n) { m_el.resize(n); }

	//! serialization
	void Serialize(DumpStream& ar);

	//! Unpack surface element data
	//! TODO: This is obsolete. Remove this.
	void UnpackLM(FEElement& el, vector<int>& lm);

public:

	//! return number of surface elements
	int Elements() { return m_el.size(); }

	//! return an element of the surface
	FESurfaceElement& Element(int i) { return m_el[i]; }

	//! returns reference to element
	FEElement& ElementRef(int n) { return m_el[n]; }

	//! find the index of a surface element
	int FindElement(FESurfaceElement& el);

public:

	//! Project a node onto a surface element
	vec3d ProjectToSurface(FESurfaceElement& el, vec3d x, double& r, double& s);

	//! check to see if a point is on element
	bool IsInsideElement(FESurfaceElement& el, double r, double s, double tol = 0);

	//! See if a ray intersects an element
	bool Intersect(FESurfaceElement& el, vec3d r, vec3d n, double rs[2], double& g, double eps);

public:
	//! calculate the surface area of a surface element
	double FaceArea(FESurfaceElement& el);

	//! return the max element size
	double MaxElementSize();

	//! calculate the metric tensor in the current configuration
	mat2d Metric(FESurfaceElement& el, double r, double s);

    //! calculate the metric tensor at an integration point
    mat2d Metric(FESurfaceElement& el, int n);
    
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

	//! calculates the covariant base vectors of a surface at an integration point
	void CoBaseVectors(FESurfaceElement& el, int j, vec3d t[2]);

	//! calculates the covariant base vectors of a surface
	void CoBaseVectors(FESurfaceElement& el, double r, double s, vec3d t[2]);

	//! calculates covariant base vectors of a surface
	void CoBaseVectors0(FESurfaceElement& el, double r, double s, vec3d t[2]);

    //! calculates contravariant base vectors of a surface  at an integration point
    void ContraBaseVectors(FESurfaceElement& el, int j, vec3d t[2]);
    
	//! calculates contravariant base vectors of a surface
	void ContraBaseVectors(FESurfaceElement& el, double r, double s, vec3d t[2]);

	//! calculates contravariant base vectors of a surface
	void ContraBaseVectors0(FESurfaceElement& el, double r, double s, vec3d t[2]);

	//! Jacobian in reference configuration for integration point n
	double jac0(FESurfaceElement& el, int n);

protected:
	vector<FESurfaceElement>	m_el;	//!< surface elements
};

#endif // !defined(AFX_FESURFACE_H__6437C4B1_5BB7_4DDA_8354_CADFF3291D3E__INCLUDED_)
