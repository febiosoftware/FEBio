#pragma once
#include "FECore/FEBoundaryCondition.h"
#include "FESurface.h"
#include "FECore/DumpFile.h"

class FESolver;

//-----------------------------------------------------------------------------
//! This is the base class for all loads that are applied to surfaces
class FESurfaceLoad : public FEBoundaryCondition
{
public:
	FESurfaceLoad(FESurface* ps);
	virtual ~FESurfaceLoad(void);

	//! Get the surface
	FESurface& Surface() { return *m_psurf; }

	//! calculate stiffness matrix
	virtual void StiffnessMatrix(FESolver* psolver) = 0;

	//! calculate residual
	virtual void Residual(FESolver* psolver, vector<double>& R) = 0;

	//! serialization
	virtual void Serialize(DumpFile& ar) = 0;

protected:
	FESurface*	m_psurf;
};
