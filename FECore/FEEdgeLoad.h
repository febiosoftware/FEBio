#pragma once
#include "FEBoundaryCondition.h"

//-----------------------------------------------------------------------------
class FEEdge;
class FEModel;
class FESolver;
class FEGlobalVector;

//-----------------------------------------------------------------------------
class FECORE_API FEEdgeLoad : public FEBoundaryCondition
{
public:
	FEEdgeLoad(FEModel* pfem);
	virtual ~FEEdgeLoad(void);

	virtual void Create(int nsegs) = 0;

	//! Set the edge to apply the load to
	void SetEdge(FEEdge* pe) { m_pedge = pe; }

	//! Get the edge
	FEEdge& Edge() { return *m_pedge; }

public:
	//! set an attribute of the surface load
	virtual bool SetAttribute(const char* szatt, const char* szval) { return false; }

public:
	//! calculate stiffness matrix
	virtual void StiffnessMatrix(FESolver* psolver) = 0;

	//! calculate residual
	virtual void Residual(FEGlobalVector& R) = 0;

protected:
	FEEdge*	m_pedge;
};
