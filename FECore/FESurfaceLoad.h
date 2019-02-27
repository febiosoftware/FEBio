#pragma once
#include "FEBoundaryCondition.h"
#include "FESurface.h"
#include "FESolver.h"
#include "FETimeInfo.h"

//-----------------------------------------------------------------------------
class FEModel;
class FEGlobalVector;

//-----------------------------------------------------------------------------
//! This is the base class for all loads that are applied to surfaces
class FECORE_API FESurfaceLoad : public FEBoundaryCondition
{
	FECORE_SUPER_CLASS

public:
	FESurfaceLoad(FEModel* pfem);
	virtual ~FESurfaceLoad(void);

	//! Set the surface to apply the load to
	virtual void SetSurface(FESurface* ps) { m_psurf = ps; }

	bool Init();

	//! Get the surface
	FESurface& GetSurface() { return *m_psurf; }

	void Serialize(DumpStream& ar) override;

public:
	//! calculate stiffness matrix
	virtual void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) {}

	//! calculate residual
	virtual void Residual(const FETimeInfo& tp, FEGlobalVector& R);

	//! unpack the surface element dofs
	virtual void UnpackLM(FESurfaceElement& el, vector<int>& lm) {}

	//! evaluate nodal values
	virtual void NodalValues(FESurfaceElement& el, vector<double>& v) {};

protected:
	FESurface*	m_psurf;
};
