#pragma once
#include "FESolver.h"
#include "DumpFile.h"
#include "FEModelComponent.h"
#include "FEGlobalVector.h"
#include "FESurface.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
// forward declaration of the model class
class FEModel;

//-----------------------------------------------------------------------------
//! Base class for nonlinear constraints enforced using an augmented Lagrangian method.

//! The constraint must provide a residual (force) contribution, its stiffness matrix
//! and an augmentation function.
//!
class FENLConstraint : public FEModelComponent
{
public:
	FENLConstraint(FEModel* pfem);
	virtual ~FENLConstraint();

public:
	virtual bool Init() = 0;
	virtual void Residual(FEGlobalVector& R) = 0;
	virtual void StiffnessMatrix(FESolver* psolver) = 0;
	virtual bool Augment(int naug) = 0;
	virtual void Serialize(DumpFile& ar) = 0;
	virtual void ShallowCopy(DumpStream& dmp, bool bsave) = 0;

	// update state
	virtual void Reset() {}
	virtual void Update() {}

	virtual FESurface* GetSurface(const char* sz) { return 0; }

	//! Get the NLC ID
	int GetID() { return m_nID; }

protected:
	int		m_nID;		//!< ID of nonlinear constraint

	static int	m_ncount;	//!< used to create unique ID's for the nonlinear constraints
};
