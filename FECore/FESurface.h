#pragma once
#include "mat2d.h"
#include "vec2d.h"
#include "FEMeshPartition.h"

//-----------------------------------------------------------------------------
class FEMesh;
class FENodeSet;
class FEFacetSet;

//-----------------------------------------------------------------------------
class FECORE_API FESurfaceMaterialPoint : public FEMaterialPoint
{
public:
	vec3d	dxr, dxs;	// tangent vectors at material point
};

//-----------------------------------------------------------------------------
//! Surface mesh

//! This class implements the basic functionality for an FE surface.
//! More specialized surfaces are derived from this class

class FECORE_API FESurface : public FEMeshPartition
{
public:
	//! default constructor
	FESurface(FEModel* fem);

	//! constructor
	FESurface(FEModel* fem, FEFacetSet* surf);

	//! destructor
	virtual ~FESurface();

	//! initialize surface data structure
	bool Init() override;
    
	//! creates surface
	void Create(int nsize, int elemType = -1) override;

	//! Build a surface from a facet set
	void BuildFromSet(FEFacetSet& set);
	void Create();

	//! serialization
	void Serialize(DumpStream& ar) override;

	//! Unpack surface element data
	//! TODO: This is obsolete. Remove this.
	void UnpackLM(FEElement& el, vector<int>& lm) override;
	
	//! Extract a node set from this surface
	// TODO: Move to MeshPartition
	FENodeSet GetNodeSet();

	//! Get a list of bools that indicate whether the corresponding node is on the boundary
	// TODO: Move to MeshPartition
	void GetBoundaryFlags(std::vector<bool>& boundary);
    
    //! Set alpha parameter for intermediate time
    void SetAlpha(const double alpha) { m_alpha = alpha; }

public:

	//! return number of surface elements
	int Elements() const override { return (int)m_el.size(); }

	//! return an element of the surface
	FESurfaceElement& Element(int i) { return m_el[i]; }

	//! return an element of the surface
	const FESurfaceElement& Element(int i) const { return m_el[i]; }

	//! returns reference to element
	FEElement& ElementRef(int n) override { return m_el[n]; }
	const FEElement& ElementRef(int n) const override { return m_el[n]; }

	//! find the solid or shell element of a surface element
	FEElement* FindElement(FESurfaceElement& el);

    //! for interface surfaces, find the index of both solid elements
    //! on either side of the interface
    void FindElements(FESurfaceElement& el);

	//! loop over all elements
	void ForEachSurfaceElement(std::function<void(FESurfaceElement& el)> f);

public:
	// Create material point data for this surface
	virtual FEMaterialPoint* CreateMaterialPoint();

	// update surface data
	void Update(const FETimeInfo& tp) override;

public:

	//! Project a node onto a surface element
	vec3d ProjectToSurface(FESurfaceElement& el, vec3d x, double& r, double& s);

	//! check to see if a point is on element
	bool IsInsideElement(FESurfaceElement& el, double r, double s, double tol = 0);

	//! See if a ray intersects an element
	bool Intersect(FESurfaceElement& el, vec3d r, vec3d n, double rs[2], double& g, double eps);

	//! Invert the surface
	void Invert();

	//! Get the spatial position given natural coordinates
	vec3d Position(FESurfaceElement& el, double r, double s);

	//! Get the nodal coordinates of an element
	void NodalCoordinates(FESurfaceElement& el, vec3d* re);

public:
	//! calculate the surface area of a surface element
	double FaceArea(FESurfaceElement& el);

	//! return the max element size
	double MaxElementSize();

	//! calculate the metric tensor in the current configuration
	mat2d Metric(FESurfaceElement& el, double r, double s);

    //! calculate the metric tensor at an integration point
    mat2d Metric(FESurfaceElement& el, int n);
    
    //! calculate the metric tensor at an integration point at previous time
    mat2d MetricP(FESurfaceElement& el, int n);
    
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

    //! calculate the global position of a point on the surface at previous time
    vec3d Local2GlobalP(FESurfaceElement& el, double r, double s);
    
    //! calculate the global position of an integration point at previous time
    vec3d Local2GlobalP(FESurfaceElement& el, int n);
    
	//! calculates the covariant base vectors of a surface at an integration point
	void CoBaseVectors(FESurfaceElement& el, int j, vec3d t[2]);

	//! calculates the covariant base vectors of a surface
	void CoBaseVectors(FESurfaceElement& el, double r, double s, vec3d t[2]);

	//! calculates covariant base vectors of a surface
	void CoBaseVectors0(FESurfaceElement& el, double r, double s, vec3d t[2]);

    //! calculates the covariant base vectors of a surface at an integration point at previoust time step
    void CoBaseVectorsP(FESurfaceElement& el, int j, vec3d t[2]);
    
    //! calculates contravariant base vectors of a surface  at an integration point
    void ContraBaseVectors(FESurfaceElement& el, int j, vec3d t[2]);
    
    //! calculates the contravariant base vectors of a surface at an integration point at previoust time step
    void ContraBaseVectorsP(FESurfaceElement& el, int j, vec3d t[2]);
    
	//! calculates contravariant base vectors of a surface
	void ContraBaseVectors(FESurfaceElement& el, double r, double s, vec3d t[2]);

	//! calculates contravariant base vectors of a surface
	void ContraBaseVectors0(FESurfaceElement& el, double r, double s, vec3d t[2]);

	//! Jacobian in reference configuration for integration point n
	double jac0(FESurfaceElement& el, int n);

	//! Jacobian in reference configuration for integration point n (and returns normal)
	double jac0(const FESurfaceElement& el, int n, vec3d& nu);

    //! Interface status
    void SetInterfaceStatus(const bool bitfc) { m_bitfc = bitfc; }
    bool GetInterfaceStatus() { return m_bitfc; }

	//! Get the facet set that created this surface
	FEFacetSet* GetFacetSet() { return m_surf; }
    
protected:
	FEFacetSet*					m_surf;		//!< the facet set from which this surface is built
	vector<FESurfaceElement>	m_el;		//!< surface elements
    bool                        m_bitfc;    //!< interface status
    double                      m_alpha;    //!< intermediate time fraction
};
