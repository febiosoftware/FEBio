#pragma once
#include "FEBoundaryCondition.h"
#include "FESurface.h"

class FESolver;

//-----------------------------------------------------------------------------
//! This is the base class for all loads that are applied to surfaces
class FESurfaceLoad : public FEBoundaryCondition
{
public:
	FESurfaceLoad(FEMesh* pm);
	virtual ~FESurfaceLoad(void);

	//! Get the surface
	FESurface& Surface() { return m_surf; }

	//! calculate stiffness matrix
	virtual void StiffnessMatrix(FESolver* psolver) = 0;

	//! calculate residual
	virtual void Residual(FESolver* psolver, vector<double>& R) = 0;

protected:
	FESurface	m_surf;
};
