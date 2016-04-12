#pragma once
#include "FEBoundaryCondition.h"
#include "FESurface.h"
#include "FESolver.h"

//-----------------------------------------------------------------------------
class FEModel;
class FEGlobalVector;

//-----------------------------------------------------------------------------
//! This is the base class for all loads that are applied to surfaces
class FESurfaceLoad : public FEBoundaryCondition
{
public:
	FESurfaceLoad(FEModel* pfem);
	virtual ~FESurfaceLoad(void);

	//! Set the surface to apply the load to
	virtual void SetSurface(FESurface* ps) { m_psurf = ps; }

	//! Get the surface
	FESurface& Surface() { return *m_psurf; }

public:
	//! set an attribute of the surface load
	virtual bool SetAttribute(const char* szatt, const char* szval) { return false; }

	//! set an attribute of a surface facet
	virtual bool SetFacetAttribute(int nface, const char* szatt, const char* szval) { return false; }

public:
	//! calculate stiffness matrix
	virtual void StiffnessMatrix(FESolver* psolver) = 0;

	//! calculate residual
	virtual void Residual(FEGlobalVector& R) = 0;

protected:
	FESurface*	m_psurf;
};
