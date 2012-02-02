#pragma once
#include "FECore/FEBoundaryCondition.h"
#include "FECore/FESurface.h"
#include "FECore/DumpFile.h"
#include "FECore/FENLSolver.h"

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
	virtual void StiffnessMatrix(FENLSolver* psolver) = 0;

	//! calculate residual
	virtual void Residual(FENLSolver* psolver, vector<double>& R) = 0;

	//! serialization
	virtual void Serialize(DumpFile& ar) = 0;

protected:
	FESurface*	m_psurf;
};
