#pragma once
#include "FEBoundaryCondition.h"
#include "FESurface.h"
#include "DumpFile.h"
#include "FESolver.h"
using namespace FECore;

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
//! This is the base class for all loads that are applied to surfaces
class FESurfaceLoad : public FEBoundaryCondition
{
public:
	FESurfaceLoad(FEModel* pfem);
	virtual ~FESurfaceLoad(void);

	virtual void Create(int nfacets) = 0;

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps) { m_psurf = ps; }

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

	//! serialization
	virtual void Serialize(DumpFile& ar) = 0;

protected:
	FESurface*	m_psurf;
	FEModel*	m_pfem;
};
